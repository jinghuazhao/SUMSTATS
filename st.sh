# 17-10-2018 JHZ

git add README.*
git commit -m "README"
for d in bmi.txt.gz st.sh
do
   git add $d
   git commit -m "$d"
done
git push
