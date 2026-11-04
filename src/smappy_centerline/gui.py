"""Accumulated ROI measurements, displayed in three linked tables."""
from PySide6.QtCore import QItemSelectionModel
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QAbstractItemView, QFileDialog, QHBoxLayout,
                              QLabel, QMessageBox, QPlainTextEdit, QPushButton,
                              QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from smappy.plugins import Result
from .plugin import display_statistics, plots


class CenterlineResultsWidget(QWidget):
    base_columns = (
        ("ROI", "ROI"), ("roi_center_x_nm", "Center X (nm)"),
        ("roi_center_y_nm", "Center Y (nm)"), ("noise_rejection", "Noise rejection"),
        ("N_before_rejection", "Input groups"),
        ("noise_excluded_nonfinite", "Invalid XYZ"),
        ("N_noise_excluded", "Noise removed"),
        ("N_after_rejection", "Retained groups"), ("noise_cutoff_nm", "Cutoff (nm)"))
    radial_columns = (
        ("ROI", "ROI"), ("radial_status", "Status"), ("dimension", "Dim."),
        ("N", "Groups"), ("centerline_length_nm", "Length (nm)"),
        ("mean_radial_distance_nm", "Mean (nm)"),
        *((f"R{p}_nm", f"R{p} (nm)") for p in (30, 50, 70, 90, 95, 100)))
    density_columns = (
        ("ROI", "ROI"), ("density_status", "Status"), ("density_dimension", "Dim."),
        ("N_grouped", "Retained groups"), ("N_grouped_in_geometry", "Groups in geometry"),
        ("geometry_measure", "Area / volume"), ("geometry_units", "Geometry units"),
        ("grouped_event_density", "Group density"), ("density_units", "Density units"),
        ("grouped_event_nyquist_nm", "Nyquist (nm)"))
    statuses = {"complete": "Complete", "disabled": "Not enabled",
                "unavailable": "Unavailable", "exploratory": "Complete (exploratory)",
                "exploratory_warning": "Model warning"}

    def __init__(self, panel):
        super().__init__(panel)
        self.panel = panel
        self.store = panel.plugin.results
        layout = QVBoxLayout(self)
        self.tables = []
        self.columns = (self.base_columns, self.radial_columns, self.density_columns)
        for title, columns in zip(("ROI information & noise rejection",
                                   "Radial Distance", "Labeling Density"), self.columns):
            label = QLabel(title)
            label.setStyleSheet("font-weight: bold")
            layout.addWidget(label)
            table = QTableWidget(0, len(columns))
            table.setHorizontalHeaderLabels([label for _, label in columns])
            table.setSelectionBehavior(QAbstractItemView.SelectRows)
            table.setSelectionMode(QAbstractItemView.ExtendedSelection)
            table.setEditTriggers(QAbstractItemView.NoEditTriggers)
            # All tables retain insertion order and the same ROI rows.
            table.setSortingEnabled(False)
            table.verticalHeader().hide()
            table.setMinimumHeight(135)
            table.itemSelectionChanged.connect(lambda t=table: self.show_selected(t))
            self.tables.append(table)
            layout.addWidget(table)
        self.table, self.radial_table, self.density_table = self.tables
        for table in self.tables:
            for other in self.tables:
                if other is not table:
                    table.verticalScrollBar().valueChanged.connect(other.verticalScrollBar().setValue)
        layout.addWidget(QLabel("Selected ROI details (all summary fields)"))
        self.details = QPlainTextEdit()
        self.details.setReadOnly(True)
        self.details.setFixedHeight(110)
        layout.addWidget(self.details)
        actions = QHBoxLayout()
        for label, callback in [("Delete selected", self.delete_selected),
                                ("Clear all", self.clear), ("Export CSV", self.export)]:
            button = QPushButton(label)
            button.clicked.connect(callback)
            actions.addWidget(button)
        layout.addLayout(actions)

    def set_result(self, result):
        self.refresh(result.data.get("centerline_roi"))

    @staticmethod
    def _text(value):
        return "—" if value is None or value == "" else (
            f"{value:.6g}" if isinstance(value, float) else str(value))

    def _values(self, entry):
        values = display_statistics(entry)
        values["radial_status"] = ("disabled" if not values.get("radial_enabled", True)
                                   else values.get("distance_status", "complete"))
        if not values.get("density_enabled", "density_status" in values):
            values["density_status"] = "disabled"
        dim = values.get("density_dimension")
        values["geometry_measure"] = values.get("Area_um2" if dim == "2D" else "Volume_um3")
        values["geometry_units"] = "µm²" if dim == "2D" else "µm³" if dim == "3D" else ""
        if values.get("noise_rejection") == "None":
            values["noise_excluded_nonfinite"] = None  # validity depends on analysis dimension
        return values

    def refresh(self, selected=None):
        for table in self.tables:
            table.blockSignals(True)
            table.clearSelection()
            table.setRowCount(len(self.store.entries))
        select = None
        for row, (identifier, entry) in enumerate(self.store.entries.items()):
            values = self._values(entry)
            for table, columns in zip(self.tables, self.columns):
                status_key = ("radial_status" if table is self.radial_table else
                              "density_status" if table is self.density_table else None)
                status = values.get(status_key)
                for col, (field, _) in enumerate(columns):
                    value = identifier if field == "ROI" else values.get(field)
                    if status == "disabled" and field not in ("ROI", status_key):
                        value = None
                    text = self.statuses.get(value, value) if field == status_key else value
                    item = QTableWidgetItem(self._text(text))
                    tooltip = self._text(value)
                    if field == status_key:
                        tooltip = values.get("distance_warning" if table is self.radial_table
                                             else "density_warnings") or self._text(text)
                        if status in ("unavailable", "exploratory_warning"):
                            item.setBackground(QColor("#ffd9d9" if status == "unavailable" else "#ffe5b4"))
                    if field == "noise_excluded_nonfinite" and values.get("noise_rejection") == "None":
                        tooltip = "No shared XYZ rejection. Invalid coordinates are counted separately by each analysis."
                    item.setToolTip(tooltip)
                    table.setItem(row, col, item)
            if identifier == selected:
                select = row
        for table in self.tables:
            table.blockSignals(False)
            table.resizeColumnsToContents()
        self.details.clear()
        if select is not None:
            self.table.selectRow(select)

    def show_selected(self, source=None):
        source = self.table if source is None else source
        rows = source.selectionModel().selectedRows()
        identifiers = {source.item(index.row(), 0).text() for index in rows}
        for table in self.tables:
            if table is source:
                continue
            table.blockSignals(True)
            table.clearSelection()
            for row in range(table.rowCount()):
                if table.item(row, 0).text() in identifiers:
                    table.selectionModel().select(table.model().index(row, 0),
                        QItemSelectionModel.Select | QItemSelectionModel.Rows)
            table.blockSignals(False)
        if not rows:
            self.details.clear()
            return
        identifier = source.item(min(index.row() for index in rows), 0).text()
        entry = self.store.entries[identifier]
        values = display_statistics(entry)
        self.details.setPlainText("\n".join(f"{key}: {self._text(value)}" for key, value in values.items()))
        self.panel.result = Result(plots=plots(entry))
        self.panel.result_button.setEnabled(bool(self.panel.result.figures()))
        if self.panel.result.figures():
            self.panel.plot()

    def delete_selected(self):
        self.store.delete([self.table.item(row.row(), 0).text()
                           for row in self.table.selectionModel().selectedRows()])
        self._after_delete()

    def clear(self):
        self.store.clear()
        self._after_delete()

    def _after_delete(self):
        self.refresh()
        self.panel.result = None
        self.panel.result_button.setEnabled(False)
        if self.panel._window is not None:
            self.panel._window.hide()

    def export(self):
        if not self.store.entries:
            self.panel.status.setText("No accumulated ROI results to export")
            return
        directory = QFileDialog.getExistingDirectory(self, "Export centerline and density CSV files")
        if not directory:
            return
        from pathlib import Path
        names = ("centerline_radial_summary.csv", "centerline_radial_all_distances.csv", "centerline_density_segments.csv",
                 "centerline_density_2d_distances.csv", "centerline_density_3d_distances.csv",
                 "centerline_summary_fields.csv")
        if any((Path(directory)/name).exists() for name in names):
            if QMessageBox.question(self, "Replace CSV files?", "Replace existing centerline CSV files?"
                                    ) != QMessageBox.Yes:
                return
        try:
            self.store.export(directory)
            self.panel.status.setText("Exported summary, distances and available density segments")
        except (OSError, ValueError) as error:
            self.panel.status.setText(str(error))


def configure_panel(panel):
    """Only the isolated extension knows about its three flat form sections."""
    for field, title in (("rejection", "1. 3D noise rejection"),
                         ("radial", "2. Radial Distance"),
                         ("density", "3. Labeling Density")):
        layout = panel.form.row_layouts[field]
        row, _ = layout.getWidgetPosition(panel.form.fields[field])
        label = QLabel(title)
        label.setStyleSheet("font-weight: bold")
        layout.insertRow(row, label)

    def update_controls(*_):
        settings = panel.form.value()
        enabled = settings.radial or settings.density
        panel.run_button.setText("Analyze/Add current ROI")
        panel.run_button.setMinimumWidth(
            panel.run_button.fontMetrics().horizontalAdvance(panel.run_button.text()) + 20)
        panel.run_button.setEnabled(enabled and not panel.busy)
        panel.run_button.setToolTip("" if enabled else "Enable at least one analysis module")
        panel.form.set_active(panel.plugin.active(settings))
        panel.form.set_shown({"maximum_distance_nm": settings.rejection == "Manual maximum distance"})

    def after_run(*_):
        if panel.busy:
            panel._thread.finished.connect(update_controls)
        else:
            update_controls()

    panel.form.changed.connect(update_controls)
    panel.ended.connect(after_run)
    panel.result_updated.connect(lambda result: panel.result_widget.set_result(result))
    update_controls()
