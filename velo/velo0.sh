samples=(
F8484/d0
F8484/d8_1
F8484/d8_2
F8492/d21_BI_1
F8492/d21_BI_2
F8492/d21_Cont_1
F8492/d21_Cont_2
F8492/d21_Nint_1
F8492/d21_Nint_2
F8655/d10_BI_1
F8655/d10_BI_2
F8655/d10_Cont_1
F8655/d10_Cont_2
F8668/d14_BI_1
F8668/d14_BI_2
F8668/d14_Cont_1
F8668/d14_Cont_2
)
ref=genes.gtf
mask=GRCm39_rmsk.gtf
for sample in ${samples[@]}; do
	sbatch -J velocyto -o log/$sample.log -p Medium -c 12 --mem 64G --wrap "velocyto run10x -m $mask -@ 12 ../../$sample/ $ref"
done
