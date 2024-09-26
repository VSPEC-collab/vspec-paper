"""
Copy stuff from static
"""
import shutil
import paths

FILENAMES = [
    'vspec-diagram.pdf'
]

for filename in FILENAMES:
    static_path = paths.static / filename
    target_path = paths.figures / filename
    shutil.copy(static_path, target_path)