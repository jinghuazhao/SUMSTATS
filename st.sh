# 18-10-2018 JHZ

git add README.*
git commit -m "README"
for d in .gitignore 97.snps bmi.txt.gz lz.sh st.sh
do
   git add $d
   git commit -m "$d"
done
git push
