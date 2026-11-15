# The names of the source files that changed since a commit, separated by a
# comma. The list feeds the variable MUTATE_FILES of the mutation script.
# Use it to run the probes of the changed code alone.
#
#     MUTATE_FILES=$(python3 test/changed_files.py HEAD~1) bash test/mutate.sh
#
# The script needs git and Python 3. It lives in the package, under test/.
# The user of the package never runs it.
import subprocess
import sys


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("usage: changed_files.py <commit> [<commit>]\n")
        return 2
    args = ["git", "diff", "--name-only"]
    args.extend(sys.argv[1:])
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    names = []
    for line in out.splitlines():
        line = line.strip()
        if not line or not line.endswith(".goose"):
            continue
        # A file under test/suite holds the tests, not the code.
        if line.startswith("test/"):
            continue
        names.append(line)
    print(",".join(names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
