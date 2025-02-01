from PyQt6.QtCore import QSize, Qt
from PyQt6.QtWidgets import QApplication, QWidget, QPushButton, QLabel, QMessageBox,  QMainWindow

# Only needed for access to command line arguments
import sys



# Subclass QMainWindow to customize your application's main window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PyTT Qt6 Interface")
        button = QPushButton("Press Me!")
        button.setCheckable(False)
        button.clicked.connect(self.quit_button_clicked)

        # Set the central widget of the Window.
        self.setCentralWidget(button)
        
    def quit_button_clicked(self):
        print("Quiting the PyTT Application")
        quit()


# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.
app = QApplication(sys.argv)

# Create a Qt widget, which will be our window.
window = MainWindow()
window.show()  # IMPORTANT!!!!! Windows are hidden by default.

# Start the event loop.
app.exec()

