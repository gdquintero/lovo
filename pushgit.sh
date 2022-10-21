now=`date +"%Y-%m-%d"`
hour=`date +"%H:%M"`
git add -A

git commit -m "Changes performed on ${now} at ${hour}"

git push
