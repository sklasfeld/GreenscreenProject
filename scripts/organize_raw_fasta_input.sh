#!/bin/bash

# create output directory
mkdir -p fastq/raw

# compress, move, and rename fastq sequence files
raw_file="SRR8180351.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR8180351.fastq > fastq/raw/inputA.fastq.gz && rm SRR8180351.fastq
fi
raw_file="SRR402844.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR402844.fastq > fastq/raw/inputB.fastq.gz && rm SRR402844.fastq
fi
raw_file="SRR6042861.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR6042861.fastq > fastq/raw/inputC.fastq.gz && rm SRR6042861.fastq
fi
raw_file="SRR8742331.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR8742331.fastq > fastq/raw/inputD.fastq.gz && rm SRR8742331.fastq
fi
raw_file="SRR5011169.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR5011169.fastq > fastq/raw/inputE.fastq.gz && rm SRR5011169.fastq
fi
raw_file="SRR8365023.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR8365023.fastq > fastq/raw/inputF.fastq.gz && rm SRR8365023.fastq
fi
raw_file="SRR6435271.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR6435271.fastq > fastq/raw/inputG.fastq.gz && rm SRR6435271.fastq
fi
raw_file="SRR7646290.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR7646290.fastq > fastq/raw/inputH.fastq.gz && rm SRR7646290.fastq
fi
raw_file="SRR8890665.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR8890665.fastq > fastq/raw/inputI.fastq.gz && rm SRR8890665.fastq
fi
raw_file1="SRR7224610.fastq"
raw_file2="SRR7224611.fastq"
raw_file3="SRR7224612.fastq"
if [ ! -f "$raw_file1"]; then
	echo "ERROR: File ${raw_file1} does not exist"
elif [ ! -f "$raw_file2"]; then
	echo "ERROR: File ${raw_file2} does not exist"
elif [ ! -f "$raw_file3"]; then
	echo "ERROR: File ${raw_file3} does not exist"
else
	cat SRR7224610.fastq SRR7224611.fastq SRR7224612.fastq > SRX4131135.fastq \
	&& rm SRR7224610.fastq SRR7224611.fastq SRR7224612.fastq
fi
raw_file="SRX4131135.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRX4131135.fastq > fastq/raw/inputJ.fastq.gz && rm SRX4131135.fastq
fi
raw_file="SRR2085420.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR2085420.fastq > fastq/raw/inputK.fastq.gz && rm SRR2085420.fastq
fi
raw_file="SRR394081.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR394081.fastq > fastq/raw/inputL.fastq.gz && rm SRR394081.fastq
fi
raw_file="SRR5681375.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR5681375.fastq > fastq/raw/inputM.fastq.gz && rm SRR5681375.fastq
fi
raw_file="SRR6768681.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR6768681.fastq > fastq/raw/inputN.fastq.gz && rm SRR6768681.fastq
fi
raw_file="SRR6989574.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR6989574.fastq > fastq/raw/inputO.fastq.gz && rm SRR6989574.fastq
fi
raw_file="SRR7264106.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR7264106.fastq > fastq/raw/inputP.fastq.gz && rm SRR7264106.fastq
fi
raw_file="SRR988541.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR988541.fastq > fastq/raw/inputQ.fastq.gz && rm SRR988541.fastq
fi
raw_file="SRR5813676.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR5813676.fastq > fastq/raw/inputR.fastq.gz && rm SRR5813676.fastq
fi
raw_file="SRR866865.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR866865.fastq > fastq/raw/inputS.fastq.gz && rm SRR866865.fastq
fi
raw_file="SRR1191642.fastq"
if [ ! -f "$raw_file"]; then
	echo "ERROR: File ${raw_file} does not exist"
else
	gzip -c SRR1191642.fastq > fastq/raw/inputT.fastq.gz && rm SRR1191642.fastq
fi
