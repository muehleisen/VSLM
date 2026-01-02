import inspect
from typing import Any
from PySide6.QtCore import QThread, Signal

class AnalysisWorker(QThread):
    """
    Executes analysis functions in a background thread to keep the GUI responsive.
    
    This worker is generic: it accepts a function object and its arguments,
    runs the function in a separate thread, and emits signals upon completion
    or error. It also supports functions that report progress via a 
    'progress_callback' argument.
    """
    
    # Emitted when the function returns successfully. Carries the return value.
    result_ready = Signal(object)
    
    # Emitted if an exception is raised. Carries the error string.
    error_occurred = Signal(str)
    
    # Emitted if the target function reports progress (0-100).
    progress_updated = Signal(int)

    def __init__(self, function: Any, *args: Any, **kwargs: Any) -> None:
        """
        Initialize the worker.
        
        Args:
            function: The python function to run.
            *args: Positional arguments to pass to the function.
            **kwargs: Keyword arguments to pass to the function.
        """
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs

    def run(self) -> None:
        """
        The main thread body. Executes the function and handles signals.
        """
        try:
            # Inspect the target function's signature.
            # If it accepts a 'progress_callback' argument, inject our emit_progress method.
            # This allows the backend logic (core.py) to update the GUI without knowing about Qt.
            sig = inspect.signature(self.function)
            if 'progress_callback' in sig.parameters:
                self.kwargs['progress_callback'] = self.emit_progress
            
            # Execute the function
            result = self.function(*self.args, **self.kwargs)
            
            # Emit the result
            self.result_ready.emit(result)
            
        except Exception as e:
            # Catch any errors and report them to the GUI
            self.error_occurred.emit(str(e))

    def emit_progress(self, value: int) -> None:
        """
        Helper method passed to the backend function to bridge the signal.
        """
        self.progress_updated.emit(value)