# bgenpy
Python module to read/write BGEN files



At the most simplist, can extract full variant data with
```python
import bgenpy
bgen = bgenpy.Reader("my.bgen")

attributes = bgen.attributes()
samples = bgen.samples()

for info in bgen:
    print(info[0])
```

## Read variant header information

If you want to just get a list of variants, then it is faster to read only the variant header information:
```python
import bgenpy
bgen = bgen.Reader("my.bgen")

bgen.seek_first_variant()
try:
    while True:
        info = bgen.read_minimal_variant()
except StopIteration:
    pass
```

## Create your own index

You can save the offsets to variants:

```python

import bgenpy
bgen = bgen.Reader("my.bgen")

bgen.seek_first_variant()
interesting_variant_offsets = []
try:
    while True:
        offset = bgen.offset()
        info = bgen.read_minimal_variant()
        if is_interesting(info):
           interesting_variant_offsets.append(offset)
except StopIteration:
    pass
```

And later jump to those variants (after re-opening the file):
```python

import bgenpy
bgen = bgen.Reader("my.bgen")

for offset in interesing_offsets:
    bgen.seek_to_variant_offset(offset)
    info = bgen.read_full_variant()
```


## zstd library

Building this module requires the zstd compression library. Then set environment variables:

```bash
export ZSTD_INC=/to/where/is/zstd-1.1.0/lib
export ZSTD_LIB=/to/where/is/zstd-1.1.0/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ZSTD_LIB

```