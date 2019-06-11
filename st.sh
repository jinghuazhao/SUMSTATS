# 11-6-2019 JHZ

git add README.*
git commit -m "README"
for d in .gitignore bmi.out 97.snps bmi.txt.gz lz.sh st.sh
do
   git add $d
   git commit -m "$d"
done
git add ST3 ST4 plasmaprotein
git commit -m "Supplementary tables"
git push
