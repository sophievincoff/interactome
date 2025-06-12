from git import Repo, InvalidGitRepositoryError
import os

def get_git_root(start_path=None):
    """
    Return the absolute path to the root of the Git repository.

    Args:
        start_path (str or None): Starting directory to search for a Git repo.
                                  Defaults to the current working directory.

    Returns:
        str or None: Path to Git root, or None if not inside a Git repo.
    """
    if start_path is None:
        start_path = os.getcwd()

    try:
        repo = Repo(start_path, search_parent_directories=True)
        return repo.git.rev_parse("--show-toplevel")
    except InvalidGitRepositoryError:
        return None

def main():
    git_root = get_git_root()
    if git_root:
        print(f"Git root: {git_root}")
    else:
        print("Not inside a Git repository.")

if __name__ == "__main__":
    main()