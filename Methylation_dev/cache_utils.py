import os
import pickle
from typing import Any, Optional, Dict
import pandas as pd
from config import PATHS, CONFIG, logger

def get_cache_filename(key: str) -> str:
    """Get cache filename with debug mode suffix if needed
    
    Args:
        key: Base name for the cache file
        
    Returns:
        Full filename including debug suffix if in debug mode
    """
    filename = f"{key}.pkl"
    if CONFIG['debug']['enabled']:
        name, ext = os.path.splitext(filename)
        return f"{name}_debug_{CONFIG['debug']['sample_size']}{ext}"
    return filename

def save_to_cache(key: str, data: Any) -> None:
    """Save data to cache file with debug mode support
    
    Args:
        key: Identifier for the cached data
        data: Data to cache (must be pickle-able)
    """
    try:
        # Create cache directory if needed
        os.makedirs(PATHS['cache_dir'], exist_ok=True)
        
        # Get cache filepath
        cache_file = get_cache_filename(key)
        cache_path = os.path.join(PATHS['cache_dir'], cache_file)
        
        # Save data
        with open(cache_path, 'wb') as f:
            pickle.dump(data, f)
        logger.info(f"Cached data saved to {cache_path}")
        
    except Exception as e:
        logger.error(f"Error saving cache for {key}: {str(e)}")

def load_from_cache(key: str) -> Optional[Any]:
    """Load data from cache file with debug mode support
    
    Args:
        key: Identifier for the cached data
        
    Returns:
        Cached data if available, None otherwise
    """
    try:
        # Get cache filepath
        cache_file = get_cache_filename(key)
        cache_path = os.path.join(PATHS['cache_dir'], cache_file)
        
        # Check if cache exists
        if not os.path.exists(cache_path):
            logger.info(f"No cache found for {key}")
            return None
            
        # Load data
        with open(cache_path, 'rb') as f:
            data = pickle.load(f)
            
        # Verify data structure for methylation results
        if key == 'methylation_results':
            if not isinstance(data, dict):
                logger.error("Cached methylation data is not a dictionary")
                return None
                
            # Verify each value is a DataFrame
            for cell_type, df in data.items():
                if not isinstance(df, pd.DataFrame):
                    logger.error(f"Cached data for {cell_type} is not a DataFrame")
                    return None
                    
        logger.info(f"Loaded data from cache: {cache_path}")
        return data
        
    except Exception as e:
        logger.error(f"Error loading cache for {key}: {str(e)}")
        return None

def clear_cache() -> None:
    """Clear all cached analysis files"""
    try:
        if os.path.exists(PATHS['cache_dir']):
            for file in os.listdir(PATHS['cache_dir']):
                if file.endswith('.pkl'):
                    os.remove(os.path.join(PATHS['cache_dir'], file))
            logger.info("Cache cleared successfully")
    except Exception as e:
        logger.error(f"Error clearing cache: {str(e)}")