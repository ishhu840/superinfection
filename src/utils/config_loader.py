"""Configuration loading utilities for the viral interference analysis system."""

import yaml
import os
from pathlib import Path
from typing import Dict, Any


class ConfigLoader:
    """Load and manage configuration settings."""
    
    def __init__(self, config_path: str = None):
        """Initialize configuration loader.
        
        Args:
            config_path: Path to configuration file. If None, uses default.
        """
        if config_path is None:
            # Default to config/config.yaml relative to project root
            project_root = Path(__file__).parent.parent.parent
            config_path = project_root / "config" / "config.yaml"
        
        self.config_path = Path(config_path)
        self._config = None
    
    def load_config(self) -> Dict[str, Any]:
        """Load configuration from YAML file.
        
        Returns:
            Dictionary containing configuration parameters.
            
        Raises:
            FileNotFoundError: If configuration file doesn't exist.
            yaml.YAMLError: If configuration file is invalid YAML.
        """
        if not self.config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {self.config_path}")
        
        try:
            with open(self.config_path, 'r') as f:
                self._config = yaml.safe_load(f)
            return self._config
        except yaml.YAMLError as e:
            raise yaml.YAMLError(f"Invalid YAML in configuration file: {e}")
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get configuration value by key.
        
        Args:
            key: Configuration key (supports dot notation, e.g., 'data_processing.min_read_length')
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        if self._config is None:
            self.load_config()
        
        keys = key.split('.')
        value = self._config
        
        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return default
    
    def get_section(self, section: str) -> Dict[str, Any]:
        """Get entire configuration section.
        
        Args:
            section: Section name (e.g., 'data_processing')
            
        Returns:
            Dictionary containing section configuration
        """
        if self._config is None:
            self.load_config()
        
        return self._config.get(section, {})
    
    def update_config(self, updates: Dict[str, Any]) -> None:
        """Update configuration with new values.
        
        Args:
            updates: Dictionary of configuration updates
        """
        if self._config is None:
            self.load_config()
        
        self._config.update(updates)
    
    def save_config(self, output_path: str = None) -> None:
        """Save current configuration to file.
        
        Args:
            output_path: Output file path. If None, overwrites original file.
        """
        if self._config is None:
            raise ValueError("No configuration loaded to save")
        
        output_path = output_path or self.config_path
        
        with open(output_path, 'w') as f:
            yaml.dump(self._config, f, default_flow_style=False, indent=2)


# Global configuration instance
_global_config = None


def get_config() -> ConfigLoader:
    """Get global configuration instance.
    
    Returns:
        Global ConfigLoader instance
    """
    global _global_config
    if _global_config is None:
        _global_config = ConfigLoader()
        _global_config.load_config()
    return _global_config


def get_config_value(key: str, default: Any = None) -> Any:
    """Convenience function to get configuration value.
    
    Args:
        key: Configuration key
        default: Default value if key not found
        
    Returns:
        Configuration value or default
    """
    return get_config().get(key, default)