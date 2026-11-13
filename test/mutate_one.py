import sys

path, old, new = sys.argv[1:4]
old = old.replace("\\&", "&").replace("\\n", "\n")
new = new.replace("\\&", "&").replace("\\n", "\n")
text = open(path).read()
if text.count(old) < 1:
    sys.exit("probe text found %d times (needs at least one) in %s: %s" % (text.count(old), path, old))
open(path, "w").write(text.replace(old, new))
