#!/bin/bash

outdir="fastq/raw"
mkdir -p ${outdir}

x="A"
raw_f="SRR8180351.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c SRR8180351.fastq > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="B"
raw_f="SRR402844.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="C"
raw_f="SRR6042861.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="D"
raw_f="SRR8742331.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="E"
raw_f="SRR5011169.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="F"
raw_f="SRR8365023.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="G"
raw_f="SRR6435271.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="H"
raw_f="SRR7646290.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="I"
raw_f="SRR8890665.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="J"
raw_f1="SRR7224610.fastq"
raw_f2="SRR7224611.fastq"
raw_f3="SRR7224612.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f1" && -f "$raw_f2" && -f "$raw_f3" && ! -f "$gzip_f" ]]; then
	cat ${raw_f1} ${raw_f2} ${raw_f3} \\
		> fastq/raw/input${x}.fastq && \
		rm ${raw_f1} ${raw_f2} ${raw_f3}
	gzip -c fastq/raw/input${x}.fastq > ${gzip_f} \\
		&& rm fastq/raw/input${x}.fastq
elif [[ ! -f "$raw_f1" || -f "$raw_f2" || -f "$raw_f3" ]]; then
	if [[ ! -f "$gzip_f"  ]]; then
		echo "ERROR: Cannot Find File(s): SRR7224610.fastq, SRR7224611.fastq, SRR7224612.fastq"
		exit
	fi
fi

x="K"
raw_f="SRR2085420.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="L"
raw_f="SRR394081.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="M"
raw_f="SRR5681375.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="N"
raw_f="SRR6768681.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="O"
raw_f="SRR6989574.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="P"
raw_f="SRR7264106.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="Q"
raw_f="SRR988541.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="R"
raw_f="SRR5813676.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="S"
raw_f="SRR866865.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi

x="T"
raw_f="SRR1191642.fastq"
gzip_f="${outdir}/input${x}.fastq.gz"
if [[ -f "$raw_f" && ! -f "$gzip_f" ]]; then
	gzip -c ${raw_f} > ${gzip_f} \\
		&& rm ${raw_f}
elif [[ ! -f "$raw_f" && ! -f "$gzip_f" ]]; then
	echo "ERROR: Cannot Find File: ${raw_f}"
	exit
fi