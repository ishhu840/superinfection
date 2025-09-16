"""Logging utilities for the viral interference analysis system."""

import logging
import logging.handlers
import os
from pathlib import Path
from typing import Optional

from .config_loader import get_config_value


def setup_logging(log_level: Optional[str] = None, 
                 log_file: Optional[str] = None,
                 max_log_size_mb: Optional[int] = None,
                 backup_count: Optional[int] = None) -> logging.Logger:
    """Set up logging configuration for the application.
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Path to log file
        max_log_size_mb: Maximum log file size in MB before rotation
        backup_count: Number of backup log files to keep
        
    Returns:
        Configured logger instance
    """
    # Get configuration values
    log_level = log_level or get_config_value('logging.level', 'INFO')
    log_file = log_file or get_config_value('logging.log_file', 'logs/viral_analysis.log')
    max_log_size_mb = max_log_size_mb or get_config_value('logging.max_log_size_mb', 10)
    backup_count = backup_count or get_config_value('logging.backup_count', 5)
    
    # Create logs directory if it doesn't exist
    log_path = Path(log_file)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Configure root logger
    logger = logging.getLogger('viral_interference')
    logger.setLevel(getattr(logging, log_level.upper()))
    
    # Remove existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler with rotation
    max_bytes = max_log_size_mb * 1024 * 1024  # Convert MB to bytes
    file_handler = logging.handlers.RotatingFileHandler(
        log_file, maxBytes=max_bytes, backupCount=backup_count
    )
    file_handler.setLevel(getattr(logging, log_level.upper()))
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    return logger


def get_logger(name: str = None) -> logging.Logger:
    """Get a logger instance.
    
    Args:
        name: Logger name. If None, returns root viral_interference logger.
        
    Returns:
        Logger instance
    """
    if name is None:
        return logging.getLogger('viral_interference')
    else:
        return logging.getLogger(f'viral_interference.{name}')


def log_function_call(func):
    """Decorator to log function calls.
    
    Args:
        func: Function to decorate
        
    Returns:
        Decorated function
    """
    def wrapper(*args, **kwargs):
        logger = get_logger(func.__module__)
        logger.debug(f"Calling {func.__name__} with args={args}, kwargs={kwargs}")
        
        try:
            result = func(*args, **kwargs)
            logger.debug(f"{func.__name__} completed successfully")
            return result
        except Exception as e:
            logger.error(f"Error in {func.__name__}: {str(e)}")
            raise
    
    return wrapper


def log_processing_step(step_name: str, logger: logging.Logger = None):
    """Context manager for logging processing steps.
    
    Args:
        step_name: Name of the processing step
        logger: Logger instance. If None, uses default logger.
        
    Example:
        with log_processing_step("Quality Control"):
            # Processing code here
            pass
    """
    class ProcessingStepLogger:
        def __init__(self, step_name: str, logger: logging.Logger):
            self.step_name = step_name
            self.logger = logger or get_logger()
        
        def __enter__(self):
            self.logger.info(f"Starting {self.step_name}...")
            return self
        
        def __exit__(self, exc_type, exc_val, exc_tb):
            if exc_type is None:
                self.logger.info(f"Completed {self.step_name} successfully")
            else:
                self.logger.error(f"Error in {self.step_name}: {exc_val}")
            return False
    
    return ProcessingStepLogger(step_name, logger)


# Initialize logging when module is imported
try:
    setup_logging()
except Exception as e:
    # Fallback to basic logging if setup fails
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logging.getLogger('viral_interference').warning(
        f"Failed to setup advanced logging: {e}. Using basic logging."
    )