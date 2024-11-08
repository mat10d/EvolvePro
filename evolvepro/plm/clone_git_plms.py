import os
import subprocess

def clone_repositories():
    # Define the base directory for external PLM repositories
    base_dir = "external"
    
    # Create the directory if it doesn't exist
    os.makedirs(base_dir, exist_ok=True)
    
    # Dictionary of repositories to clone
    repositories = {
        "proteinbert": {
            "url": "https://github.com/nadavbra/protein_bert.git",
            "post_clone": [
                ["git", "submodule", "init"],
                ["git", "submodule", "update"],
                ["python", "setup.py", "install"]
            ]
        },
        "efficient-evolution": {
            "url": "https://github.com/brianhie/efficient-evolution.git",
            "post_clone": []  # No post-clone commands needed
        }
    }
    
    # Clone each repository and run post-clone commands
    for folder_name, repo_info in repositories.items():
        target_dir = os.path.join(base_dir, folder_name)
        if not os.path.exists(target_dir):
            print(f"Cloning {repo_info['url']} into {target_dir}...")
            subprocess.run(["git", "clone", repo_info['url'], target_dir])
            
            # Change to the target directory for post-clone commands
            original_dir = os.getcwd()
            os.chdir(target_dir)
            
            # Run any post-clone commands
            for cmd in repo_info['post_clone']:
                print(f"Running: {' '.join(cmd)}")
                subprocess.run(cmd)
            
            # Change back to original directory
            os.chdir(original_dir)
        else:
            print(f"Directory {target_dir} already exists, skipping...")

if __name__ == "__main__":
    clone_repositories()