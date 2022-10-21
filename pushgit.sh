now=`date +"%m-%d-%Y"`
hour=`date +"%H:%M"`
git add -A

git commit -m "Changes performed on ${now} at ${hour}"

git push
