import os
import git
import hashlib
import datetime
import json
import shutil
import glob
from pathlib import Path
import tempfile
from PIL import Image
import sys

def extract_git_images(repo_path, output_dir):
    """
    Extract images from all branches and commits in a local Git repository.

    Args:
        repo_path: Path to the local Git repository
        output_dir: Directory to store the extracted images
    """
    print(f"Starting image extraction from local repository at {repo_path}")

    # Validate repo path
    repo_path_obj = Path(repo_path)
    if not (repo_path_obj / '.git').is_dir():
        print(f"Error: '{repo_path}' is not a valid Git repository.")
        return None

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    try:
        repo = git.Repo(repo_path)
    except git.InvalidGitRepositoryError:
        print(f"Error: Invalid Git repository at {repo_path}")
        return None
    except Exception as e:
        print(f"Error opening repository: {e}")
        return None

    # Create results directory structure
    results_dir = output_path / 'images'
    results_dir.mkdir(exist_ok=True)

    # Create metadata directory
    metadata_dir = output_path / 'metadata'
    metadata_dir.mkdir(exist_ok=True)

    # Get all local branches
    branches = [head.name for head in repo.heads]
    # Also include remote branches if needed, requires fetching first
    try:
        repo.git.fetch('--all')
        remote_branches = []
        for ref in repo.remote().refs:
            if 'HEAD' not in ref.name:
                 # Format remote branch name (e.g., origin/main -> main)
                 # Handle cases where branch names might contain '/'
                 branch_name_parts = ref.name.split('/')
                 if len(branch_name_parts) > 1:
                     local_branch_name = '/'.join(branch_name_parts[1:])
                     if local_branch_name not in branches:
                         remote_branches.append(ref.name) # Store the full remote ref name
                         branches.append(local_branch_name) # Add the simplified name for processing if not already present locally
                 elif ref.name not in branches: # Handle cases like just 'main' from remote? Unlikely but possible.
                     branches.append(ref.name)

        # Add unique remote branches to the list to iterate over
        # We need the full ref name (e.g., 'origin/main') to get commits correctly
        all_refs_to_check = [head.name for head in repo.heads] + remote_branches
        all_refs_to_check = list(set(all_refs_to_check)) # Ensure uniqueness

    except git.GitCommandError as e:
        print(f"Warning: Could not fetch remote branches: {e}. Processing local branches only.")
        all_refs_to_check = [head.name for head in repo.heads]


    print(f"Found {len(all_refs_to_check)} refs (local branches and remote refs) to check.")

    # Dictionary to store unique images
    unique_images = {}

    # Track overall statistics
    stats = {
        'total_refs': len(all_refs_to_check),
        'processed_refs': 0,
        'total_commits': 0,
        'total_images': 0,
        'unique_images': 0
    }

    # Store current branch to return later
    try:
        original_branch = repo.active_branch.name
    except TypeError: # Detached HEAD state
        original_branch = repo.head.commit.hexsha
        print("Warning: Repository is in a detached HEAD state.")


    # Process each ref (branch or remote ref)
    processed_commits = set() # Avoid processing the same commit multiple times if it's in multiple branches

    for ref_name in all_refs_to_check:
        local_branch_name = ref_name.split('/')[-1] # Simplified name for directory structure
        print(f"\nProcessing ref: {ref_name} (as {local_branch_name})")
        stats['processed_refs'] += 1

        # Create a directory for this branch/ref
        # Replace characters invalid for directory names
        safe_branch_name = local_branch_name.replace('/', '_').replace('\\', '_')
        branch_dir = results_dir / safe_branch_name
        branch_dir.mkdir(exist_ok=True)

        # Get commit history for this ref
        try:
            commits = list(repo.iter_commits(ref_name))
            print(f"  Found {len(commits)} commits for this ref")
        except git.GitCommandError as e:
            print(f"  Error getting commit history for {ref_name}: {e}")
            continue
        except Exception as e:
            print(f"  Unexpected error getting commits for {ref_name}: {e}")
            continue


        # Process each commit
        for commit in commits:
            commit_hash = commit.hexsha
            if commit_hash in processed_commits:
                # print(f"  Skipping already processed commit: {commit_hash[:8]}")
                continue

            processed_commits.add(commit_hash)
            stats['total_commits'] += 1 # Count unique commits processed

            short_hash = commit_hash[:8]
            commit_date = datetime.datetime.fromtimestamp(commit.committed_date)
            commit_date_str = commit_date.strftime('%Y%m%d_%H%M%S')

            print(f"  Processing commit: {short_hash} from {commit_date_str}")

            # Create directory for this commit
            commit_dir = branch_dir / f"{commit_date_str}_{short_hash}"

            # --- Efficiently access files in the commit ---
            try:
                tree = commit.tree
            except Exception as e:
                print(f"    Error accessing tree for commit {short_hash}: {e}")
                continue

            image_extensions = ['.png', '.jpg', '.jpeg', '.gif', '.bmp', '.tiff', '.webp', '.svg']
            commit_images_found = 0

            # Recursively traverse the tree
            for blob in tree.traverse():
                if blob.type == 'blob' and Path(blob.path).suffix.lower() in image_extensions:
                    rel_path = Path(blob.path)
                    # print(f"    Found image: {rel_path}") # Debugging
                    commit_images_found += 1
                    stats['total_images'] += 1 # Count every instance found

                    # Create the commit directory only if we find images
                    commit_dir.mkdir(exist_ok=True)

                    # Get blob content and calculate hash first
                    try:
                        content = blob.data_stream.read()
                        file_hash = hashlib.md5(content).hexdigest()
                    except Exception as e:
                        print(f"    Error reading blob data for {rel_path} in commit {short_hash}: {e}")
                        continue

                    # Image destination path (flattened using hash)
                    dest_file = commit_dir / f"{file_hash}{rel_path.suffix}"
                    # commit_dir is already created if images are found

                    # Store image metadata if it's a new unique image
                    if file_hash not in unique_images:
                        dimensions = None
                        # Try to get image dimensions using Pillow from memory
                        try:
                            from io import BytesIO
                            with Image.open(BytesIO(content)) as img:
                                dimensions = img.size
                        except Exception as img_err:
                            # print(f"      Could not get dimensions for {rel_path}: {img_err}") # Debugging
                            pass # Keep dimensions as None if PIL fails

                        unique_images[file_hash] = {
                            'hash': file_hash,
                            'original_path': str(rel_path), # Store the first path found
                            'size_bytes': len(content),
                            'dimensions': dimensions,
                            'occurrences': []
                        }
                        stats['unique_images'] += 1

                    # Add occurrence details
                    occurrence_info = {
                        'branch': local_branch_name, # Use simplified name
                        'ref': ref_name, # Store the full ref name checked
                        'commit': commit_hash,
                        'commit_short': short_hash,
                        'date': commit_date.isoformat(),
                        'path': str(rel_path) # Path within this specific commit
                    }

                    # Avoid duplicate occurrences for the same image in the same commit/branch combo
                    # This can happen if the same commit is reachable from multiple refs we check
                    existing_occurrences = unique_images[file_hash]['occurrences']
                    is_duplicate_occurrence = any(
                        occ['commit'] == commit_hash and occ['path'] == str(rel_path) and occ['ref'] == ref_name
                        for occ in existing_occurrences
                    )

                    if not is_duplicate_occurrence:
                         unique_images[file_hash]['occurrences'].append(occurrence_info)


                    # Write the blob content to the destination file
                    try:
                        with open(dest_file, 'wb') as f_out:
                            f_out.write(content)
                    except Exception as e:
                        print(f"    Error writing image file {dest_file}: {e}")

            if commit_images_found > 0:
                 print(f"    Found {commit_images_found} images in commit {short_hash}")
            # else:
            #      print(f"    No images found in commit {short_hash}") # Debugging

    # --- Restore original state ---
    print("\nRestoring original repository state...")
    try:
        repo.git.checkout(original_branch)
        print(f"Checked out original ref: {original_branch}")
    except git.GitCommandError as e:
        print(f"Error checking out original ref '{original_branch}': {e}")
        print("Repository might be in an unclean state.")
    except Exception as e:
         print(f"Unexpected error restoring state: {e}")


    # Save metadata
    metadata_file = metadata_dir / 'images.json'
    print(f"\nSaving metadata to {metadata_file}...")
    try:
        with open(metadata_file, 'w') as f:
            json.dump(unique_images, f, indent=2)
    except Exception as e:
        print(f"Error saving metadata: {e}")

    # Create an index file with search capabilities
    html_index_file = output_path / 'index.html'
    print(f"Creating HTML index at {html_index_file}...")
    try:
        create_html_index(output_path, unique_images, stats)
    except Exception as e:
        print(f"Error creating HTML index: {e}")

    print(f"\nExtraction complete!")
    print(f"Processed {stats['processed_refs']} refs and {stats['total_commits']} unique commits")
    print(f"Found {stats['unique_images']} unique images out of {stats['total_images']} total image instances found across history")
    print(f"Results stored in: {output_path.resolve()}")

    return output_dir

def create_html_index(output_path, unique_images, stats):
    """Create an HTML index page with search functionality"""
    results_dir = output_path / 'images'
    html_file_path = output_path / 'index.html'

    with open(html_file_path, 'w', encoding='utf-8') as f:
        f.write('''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Git Repository Image Search</title>
    <style>
        body { font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell, "Open Sans", "Helvetica Neue", sans-serif; margin: 0; padding: 20px; background-color: #f8f9fa; color: #212529; }
        .container { max-width: 1600px; margin: 0 auto; }
        h1, h2 { color: #343a40; }
        .stats { background-color: #e9ecef; padding: 15px; border-radius: 8px; margin-bottom: 20px; border: 1px solid #dee2e6; }
        .controls { background-color: #fff; padding: 20px; border-radius: 8px; margin-bottom: 20px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .search, .filters { margin-bottom: 15px; }
        label { display: block; margin-bottom: 5px; font-weight: 500; color: #495057; }
        input[type="text"], select { width: 100%; padding: 10px; box-sizing: border-box; font-size: 14px; border: 1px solid #ced4da; border-radius: 4px; }
        input[type="text"]:focus, select:focus { border-color: #80bdff; outline: 0; box-shadow: 0 0 0 0.2rem rgba(0,123,255,.25); }
        .filter-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; }
        .image-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 20px; }
        .image-item { background-color: #fff; border: 1px solid #dee2e6; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.05); transition: box-shadow 0.2s ease-in-out; }
        .image-item:hover { box-shadow: 0 4px 8px rgba(0,0,0,0.1); }
        .image-item img { max-width: 100%; display: block; height: 200px; object-fit: contain; margin: 0 auto; background: #f8f9fa; border-bottom: 1px solid #dee2e6; }
        .image-info { padding: 15px; font-size: 13px; line-height: 1.5; }
        .image-info p { margin: 0 0 8px 0; word-wrap: break-word; }
        .image-info strong { color: #495057; }
        .image-info .path { font-family: monospace; font-size: 12px; color: #6c757d; }
        .image-info a { color: #007bff; text-decoration: none; }
        .image-info a:hover { text-decoration: underline; }
        .occurrences-toggle { cursor: pointer; color: #007bff; font-weight: bold; margin-top: 5px; display: inline-block; }
        .occurrences-list { display: none; margin-top: 8px; padding-left: 15px; border-left: 2px solid #e9ecef; font-size: 12px; max-height: 150px; overflow-y: auto; }
        .occurrences-list li { margin-bottom: 5px; }
        .no-results { text-align: center; padding: 40px; color: #6c757d; font-size: 1.2em; }
    </style>
</head>
<body>
    <div class="container">
        <h1>Git Repository Image Search</h1>

        <div class="stats">
            <h2>Statistics</h2>
            <p>Refs processed: ${stats['processed_refs']} / ${stats['total_refs']}</p>
            <p>Unique commits processed: ${stats['total_commits']}</p>
            <p>Total image instances found: ${stats['total_images']}</p>
            <p>Unique images found: ${stats['unique_images']}</p>
        </div>

        <div class="controls">
            <div class="search">
                <label for="searchInput">Search by filename or path:</label>
                <input type="text" id="searchInput" placeholder="e.g., logo.png or assets/images">
            </div>
            <div class="filter-grid">
                <div class="filter-group">
                    <label for="branchFilter">Branch/Ref:</label>
                    <select id="branchFilter">
                        <option value="all">All Branches/Refs</option>
                    </select>
                </div>
                <div class="filter-group">
                    <label for="sortBy">Sort by:</label>
                    <select id="sortBy">
                        <option value="path">File Path (A-Z)</option>
                        <option value="path-desc">File Path (Z-A)</option>
                        <option value="date-desc">Last Seen (Newest First)</option>
                        <option value="date-asc">First Seen (Oldest First)</option>
                        <option value="occurrences-desc">Occurrences (Most First)</option>
                        <option value="occurrences-asc">Occurrences (Fewest First)</option>
                        <option value="size-desc">Size (Largest First)</option>
                        <option value="size-asc">Size (Smallest First)</option>
                    </select>
                </div>
            </div>
        </div>

        <div class="image-grid" id="imageGrid">
''')

        # Collect all unique branches/refs from occurrences for the filter dropdown
        all_branches_refs = set()
        image_data_for_js = []

        # Prepare data for HTML and JavaScript
        for file_hash, metadata in unique_images.items():
            if not metadata['occurrences']: continue # Skip if no occurrences somehow

            # Sort occurrences by date, newest first
            occurrences = sorted(metadata['occurrences'], key=lambda x: x['date'], reverse=True)
            latest_occurrence = occurrences[0]
            earliest_occurrence = occurrences[-1]

            # Use the path from the latest occurrence for display image
            latest_path = Path(latest_occurrence['path'])
            latest_branch_dir_name = latest_occurrence['branch'].replace('/', '_').replace('\\', '_')
            latest_commit_str = f"{datetime.datetime.fromisoformat(latest_occurrence['date']).strftime('%Y%m%d_%H%M%S')}_{latest_occurrence['commit_short']}"

            # Construct the relative path to the image copy (using hash filename)
            # The actual saved file is named by hash + suffix
            hashed_filename = f"{file_hash}{latest_path.suffix}"
            relative_image_path = Path('images') / latest_branch_dir_name / latest_commit_str / hashed_filename

            # Ensure POSIX-style paths for HTML src attribute
            display_image_src = relative_image_path.as_posix()

            dimensions_str = f"{metadata['dimensions'][0]}x{metadata['dimensions'][1]}" if metadata['dimensions'] else "N/A"
            size_kb = metadata['size_bytes'] / 1024

            # Collect unique branches/refs for the filter
            for occ in occurrences:
                all_branches_refs.add(occ['branch']) # Use the simplified branch name for filtering

            # Data for sorting/filtering in JS
            image_item_data = {
                'hash': file_hash,
                'path': metadata['original_path'], # Use the canonical original path
                'branches': list(set(occ['branch'] for occ in occurrences)), # Simplified names
                'occurrences_count': len(occurrences),
                'first_seen': earliest_occurrence['date'],
                'last_seen': latest_occurrence['date'],
                'size_bytes': metadata['size_bytes']
            }
            image_data_for_js.append(image_item_data)

            # Write HTML for the image item
            f.write(f'''
        <div class="image-item"
             data-hash="{file_hash}"
             data-path="{metadata['original_path'].lower()}"
             data-branches='{json.dumps(list(set(occ['branch'] for occ in occurrences)))}'
             data-occurrences="{len(occurrences)}"
             data-first-seen="{earliest_occurrence['date']}"
             data-last-seen="{latest_occurrence['date']}"
             data-size="{metadata['size_bytes']}">
            <a href="{display_image_src}" target="_blank" title="Click to open latest version in new tab">
                <img src="{display_image_src}" alt="{latest_path.name}" loading="lazy">
            </a>
            <div class="image-info">
                <p><strong>File:</strong> {latest_path.name}</p>
                <p><strong>Path:</strong> <span class="path">{latest_path.parent}</span></p>
                <p><strong>Size:</strong> {size_kb:.1f} KB</p>
                <p><strong>Dimensions:</strong> {dimensions_str}</p>
                <p><strong>Unique Occurrences:</strong> {len(occurrences)}</p>
                <p><strong>First Seen:</strong> {datetime.datetime.fromisoformat(earliest_occurrence['date']).strftime('%Y-%m-%d')}</p>
                <p><strong>Last Seen:</strong> {datetime.datetime.fromisoformat(latest_occurrence['date']).strftime('%Y-%m-%d')}</p>
                <details>
                    <summary class="occurrences-toggle">Show History ({len(occurrences)})</summary>
                    <ul class="occurrences-list">''')
            # Add list of occurrences
            for occ in occurrences:
                 occ_date_str = datetime.datetime.fromisoformat(occ['date']).strftime('%Y-%m-%d %H:%M')
                 # Construct link to the specific version of the image
                 occ_branch_dir = occ['branch'].replace('/', '_').replace('\\', '_')
                 occ_commit_str = f"{datetime.datetime.fromisoformat(occ['date']).strftime('%Y%m%d_%H%M%S')}_{occ['commit_short']}"
                 # Link to the hashed filename for this specific occurrence
                 occ_hashed_filename = f"{file_hash}{Path(occ['path']).suffix}"
                 occ_img_path = Path('images') / occ_branch_dir / occ_commit_str / occ_hashed_filename
                 f.write(f'<li><a href="{occ_img_path.as_posix()}" target="_blank">{occ_date_str}</a> ({occ["commit_short"]})<br>Ref: {occ["ref"]}<br>Path: {occ["path"]}</li>')

            f.write('''
                    </ul>
                </details>
            </div>
        </div>''')

        # Close the image grid and add JavaScript
        f.write('''
        </div>
        <div id="noResults" class="no-results" style="display: none;">No images match your criteria.</div>
    </div>

    <script>
        const imageGrid = document.getElementById('imageGrid');
        const imageItems = Array.from(imageGrid.querySelectorAll('.image-item'));
        const searchInput = document.getElementById('searchInput');
        const branchFilter = document.getElementById('branchFilter');
        const sortBySelect = document.getElementById('sortBy');
        const noResultsDiv = document.getElementById('noResults');

        // Populate branch filter
        const branches = ''')
        f.write(json.dumps(sorted(list(all_branches_refs))))
        f.write(''';
        branches.forEach(branch => {
            const option = document.createElement('option');
            option.value = branch;
            option.textContent = branch;
            branchFilter.appendChild(option);
        });

        function updateDisplay() {
            const searchValue = searchInput.value.toLowerCase().trim();
            const branchValue = branchFilter.value;
            const sortValue = sortBySelect.value;

            // Filter images
            let visibleImages = imageItems.filter(item => {
                const itemPath = item.dataset.path;
                const itemBranches = JSON.parse(item.dataset.branches); // Branches are stored as JSON string

                const matchesSearch = !searchValue || itemPath.includes(searchValue);
                const matchesBranch = branchValue === 'all' || itemBranches.includes(branchValue);

                // Hide or show item based on filter
                const isVisible = matchesSearch && matchesBranch;
                item.style.display = isVisible ? 'block' : 'none';
                return isVisible;
            });

            // Sort visible images
            visibleImages.sort((a, b) => {
                const pathA = a.dataset.path;
                const pathB = b.dataset.path;
                const lastSeenA = new Date(a.dataset.lastSeen);
                const lastSeenB = new Date(b.dataset.lastSeen);
                const firstSeenA = new Date(a.dataset.firstSeen);
                const firstSeenB = new Date(b.dataset.firstSeen);
                const occurrencesA = parseInt(a.dataset.occurrences);
                const occurrencesB = parseInt(b.dataset.occurrences);
                const sizeA = parseInt(a.dataset.size);
                const sizeB = parseInt(b.dataset.size);

                switch(sortValue) {
                    case 'path': return pathA.localeCompare(pathB);
                    case 'path-desc': return pathB.localeCompare(pathA);
                    case 'date-desc': return lastSeenB - lastSeenA; // Newest first
                    case 'date-asc': return firstSeenA - firstSeenB; // Oldest first
                    case 'occurrences-desc': return occurrencesB - occurrencesA;
                    case 'occurrences-asc': return occurrencesA - occurrencesB;
                    case 'size-desc': return sizeB - sizeA;
                    case 'size-asc': return sizeA - sizeB;
                    default: return 0;
                }
            });

            // Reorder DOM elements
            visibleImages.forEach(item => {
                imageGrid.appendChild(item); // Move item to the end of the grid container
            });

            // Show/hide 'no results' message
            noResultsDiv.style.display = visibleImages.length === 0 ? 'block' : 'none';
        }

        // Debounce function
        function debounce(func, wait) {
            let timeout;
            return function executedFunction(...args) {
                const later = () => {
                    clearTimeout(timeout);
                    func(...args);
                };
                clearTimeout(timeout);
                timeout = setTimeout(later, wait);
            };
        }

        // Event listeners
        searchInput.addEventListener('input', debounce(updateDisplay, 300)); // Debounce search input
        branchFilter.addEventListener('change', updateDisplay);
        sortBySelect.addEventListener('change', updateDisplay);

        // Initial display
        updateDisplay();
    </script>
</body>
</html>''')

if __name__ == "__main__":
    import argparse

    # Use the current directory if no repo path is provided
    default_repo_path = '.'
    # Suggest an output directory within the script's location or CWD
    default_output_dir = './git_image_results'

    parser = argparse.ArgumentParser(description='Extract images from a Git repository history.')
    parser.add_argument('repo_path', nargs='?', default=default_repo_path,
                        help=f'Path to the local Git repository (default: {default_repo_path})')
    parser.add_argument('--output', default=default_output_dir,
                        help=f'Directory to store results (default: {default_output_dir})')

    args = parser.parse_args()

    # Use absolute paths for clarity
    repo_abs_path = Path(args.repo_path).resolve()
    output_abs_path = Path(args.output).resolve()


    print(f"Repository to analyze: {repo_abs_path}")
    print(f"Output directory: {output_abs_path}")

    extract_git_images(str(repo_abs_path), str(output_abs_path))
