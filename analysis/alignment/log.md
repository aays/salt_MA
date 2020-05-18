
## 18/5/2020

today:

- downloading all these sweet, sweet fastq files

command in private notion project page, for obvious reasons - needed to use
`wget --no-remove-listing` to get the actual file listings and then `wget -r` to
ensure recursive downloading

took 30 minutes but it's looks like we're good to go! 

checking checksums:

```python
import subprocess
from tqdm import tqdm
import glob
import os
import re

counter = 0
pass_counter = 0
dirnames = os.listdir()
for d in tqdm(dirnames):
    if os.path.isdir(d) and re.search('[A-Z]{2}[0-9]{2,4}_[05]', d):
        fastq_files = glob.glob(d + '/*fq.gz')
        checksum_file = d + '/MD5.txt'
        with open(checksum_file, 'r') as f:
            lines = [line.rstrip('\n').split('  ') for line in f] # two spaces
            checksums = {d + '/' + line[1]: line[0] for line in lines} # fname: checksum
        own_checksums = dict.fromkeys(fastq_files)
        for fname in fastq_files:
            proc = subprocess.Popen(['md5sum', fname],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            ch_split = out.decode('utf-8').rstrip('\n').split('  ')
            if checksums[fname] == ch_split[0]:
                pass_counter += 1
            counter += 1

print(counter, pass_counter)
# 82 82
```

finally, changing permissions for `data/fastq` to 544 just to be safe (dirs need
to be set to r-x to even be openable)
