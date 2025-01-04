import os
import pickle
from typing import Any, Optional
from config import PATHS, CONFIG, logger

def get_cache_filename(base_filename: str) -> str:
    """Get cache filename with debug mode suffix if needed"""
    if CONFIG['debug']['enabled']:
        name, ext = os.path.splitext(base_filename)
        return f"{name}_debug_{CONFIG['debug']['sample_size']}{ext}"
    return base_filename

def save_to_cache(data: Any, cache_file: str) -> None:
    """Save data to cache file with debug mode support"""
    try:
        os.makedirs(PATHS['cache_dir'], exist_ok=True)
        cache_file = get_cache_filename(cache_file)
        cache_path = os.path.join(PATHS['cache_dir'], cache_file)
        with open(cache_path, 'wb') as f:
            pickle.dump(data, f)
        logger.info(f"Cached data saved to {cache_path}")
    except Exception as e:
        logger.error(f"Error saving cache: {str(e)}")

def load_from_cache(cache_file: str) -> Optional[Any]:
    """Load data from cache file with debug mode support"""
    try:
        cache_file = get_cache_filename(cache_file)
        cache_path = os.path.join(PATHS['cache_dir'], cache_file)
        if os.path.exists(cache_path):
            with open(cache_path, 'rb') as f:
                data = pickle.load(f)
            logger.info(f"Loaded data from cache: {cache_path}")
            return data
    except Exception as e:
        logger.error(f"Error loading cache: {str(e)}")
    return None

def clear_cache() -> None:
    """Clear all cached analysis files"""
    try:
        if os.path.exists(PATHS['cache_dir']):
            for file in os.listdir(PATHS['cache_dir']):
                os.remove(os.path.join(PATHS['cache_dir'], file))
            logger.info("Cache cleared successfully")
    except Exception as e:
        logger.error(f"Error clearing cache: {str(e)}")