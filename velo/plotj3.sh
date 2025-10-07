for name in d{10,14,21}_{BI,Cont}; do
	python3 plotj3.py --params params.toml --name $name result/${name}_1/${name}_1.h5ad result/${name}_2/${name}_2.h5ad
done
