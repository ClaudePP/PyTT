############################################################
#                    PyTT MAIN PROGRAM                     #
############################################################

# --------- Importing the necessary Libraries  ----------- #

import sys
import os
#from Modules import TargetGeometry 
#from Modules import FileIO
#from Modules import CoreSimulationFunctions
#from Modules import NecessaryVariables as nv



import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
)
        
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget,
    QPushButton,
    QLineEdit,
    QLabel,
    QHBoxLayout,
    QMessageBox,
)

class PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig, self.ax = plt.subplots()
        super().__init__(self.fig)
        self.setParent(parent)

    def plot(self, x, y):
        self.ax.clear()
        self.ax.plot(x, y, label="Line Plot")
        self.ax.set_title("Sample Plot")
        self.ax.set_xlabel("X-axis")
        self.ax.set_ylabel("Y-axis")
        self.ax.legend()
        self.draw()

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PyQt6 Plot Example")
        self.setGeometry(100, 100, 800, 600)

        # Central widget and layout
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)

        # Input section
        self.input_layout = QHBoxLayout()
        self.layout.addLayout(self.input_layout)

        self.input_x = QLineEdit()
        self.input_x.setPlaceholderText("Enter X values (comma-separated)")
        self.input_layout.addWidget(QLabel("X Values:"))
        self.input_layout.addWidget(self.input_x)

        self.input_y = QLineEdit()
        self.input_y.setPlaceholderText("Enter Y values (comma-separated)")
        self.input_layout.addWidget(QLabel("Y Values:"))
        self.input_layout.addWidget(self.input_y)

        # Plot button
        self.plot_button = QPushButton("Generate Plot")
        self.plot_button.clicked.connect(self.plot_graph)
        self.layout.addWidget(self.plot_button)

        # Plot canvas
        self.canvas = PlotCanvas(self)
        self.layout.addWidget(self.canvas)

    def plot_graph(self):
        try:
            # Parse input values
            x = list(map(float, self.input_x.text().split(",")))
            y = list(map(float, self.input_y.text().split(",")))

            if len(x) != len(y):
                raise ValueError("X and Y values must have the same length!")

            # Plot the graph
            self.canvas.plot(x, y)
        except ValueError as e:
            QMessageBox.critical(self, "Error", str(e))

        
# The code starts here, when the user executes this MAIN.py file. 
# Run the application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())

