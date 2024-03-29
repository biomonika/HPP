#use palindrover results to generate palindrome track

rm chrY.palindromes.bed
sed 's/chrY/hs1/' hg002.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs2/' mPanTro3.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs3/' mPanPan1.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs4/' mGorGor1.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs5/' mPonAbe1.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs6/' mPonPyg2.chrY.palindromes.bed >>chrY.palindromes.bed
sed 's/chrY/hs7/' mSymSyn1.chrY.palindromes.bed >>chrY.palindromes.bed

for a in *palindromes.bed; do echo $a; python merge_palindrover_results.py $a; done;
bedtools merge -i chrY.palindromes.bed.palindrover.bed >palindrover.merged.bed
