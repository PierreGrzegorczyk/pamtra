time for f in $(git ls-files); do
    echo "Adding $f"
    /usr/bin/time -f "%e s" git add "$f"
done
