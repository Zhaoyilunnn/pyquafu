
"""
Exceptions for errors raised while building circuit.
"""

class IqcsError(Exception):
    """Base class for errors raised by Iqcs."""

    def __init__(self, *message):
        """Set the error message."""
        super().__init__(" ".join(message))
        self.message = " ".join(message)

    def __str__(self):
        """Return the message."""
        return repr(self.message)

class CircuitError(IqcsError):
    """Exceptions for errors raised while building circuit."""
    pass

class ServerError(IqcsError):
    pass

class CompileError(IqcsError):
    pass

class UserError(IqcsError):
    pass