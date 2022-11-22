import os
import re
import tempfile
import subprocess
from pathlib import Path
import tarfile
import gzip
import shutil

## VIASH START
par = {
  'genomeDir': 'resources_test/cellranger_tiny_fastq/cellranger_tiny_ref/star/',
  'readFilesIn': [
    'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R1_001.fastq.gz',
    'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L001_R2_001.fastq.gz',
    'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R1_001.fastq.gz',
    'resources_test/cellranger_tiny_fastq/cellranger_tiny_fastq/tinygex_S1_L002_R2_001.fastq.gz'
  ]
}
meta = {
  'config': 'src/mapping/star_align/.config.vsh.yaml',
  'cpus': None,
  'memory_b': None,
  'memory_kb': None,
  'memory_mb': None,
  'memory_gb': None,
  'memory_tb': None,
  'memory_pb': None,
  'temp_dir': '/tmp'
}
## VIASH END

# fastqgz_regex = r'([A-Za-z0-9\-_\.]+)_S(\d+)_L(\d+)_([RI]\d+)_(\d+)\.fastq(\.gz)?'
fastqgz_regex = r'([A-Za-z0-9\-_\.]+)_S(\d+)_L(\d+)_(R\d+)_(\d+)\.fastq(\.gz)?'

def is_gz_file(filepath):
  with open(filepath, 'rb') as test_f:
    return test_f.read(2) == b'\x1f\x8b'

# rename keys
to_rename = {'input': 'readFilesIn', 'reference': 'genomeDir', 'output': 'outFileNamePrefix'}
par = { to_rename[key] if key in to_rename.keys() else key: value
  for key, value in par.items() }

# create output dir if need be
if not os.path.exists(par["outFileNamePrefix"]):
  os.mkdir(par["outFileNamePrefix"])

# make sure there is a trailing /
par["outFileNamePrefix"] = par["outFileNamePrefix"] + "/"

with tempfile.TemporaryDirectory(prefix="star-", dir=meta["temp_dir"]) as temp_dir:
  print(">> Check whether input files are directories", flush=True)
  new_read_files_in = []
  for value in par["readFilesIn"]:
    value_path = Path(value)
    if value_path.is_dir():
      print(f"Input '{value_path}' is a directory, traversing to see if we can detect any FASTQ files.", flush=True)
      value_paths = [file for file in value_path.rglob('*') if re.match(fastqgz_regex, file.name) ]
      new_read_files_in.extend(value_paths)
    else:
      new_read_files_in.append(value_path)
    par["readFilesIn"] = new_read_files_in
  print("", flush=True)

  # checking for compressed files, ungzip files if need be
  temp_dir_path = Path(temp_dir)
  input_pars = ["genomeDir", "readFilesIn"]
  for par_name in input_pars:
    par_values = par[par_name]

    # turn value into list
    is_multiple = isinstance(par_values, list)
    if not is_multiple:
      par_values = [ par_values ]
    
    # output list
    new_values = []
    
    for par_value in par_values:
      print(f'>> Check compression of --{par_name} with value: {par_value}', flush=True)
      if Path(par_value).is_file() and tarfile.is_tarfile(par_value):
        # Remove two extensions (if they exist)
        extaction_dir_name = Path(par_value.stem).stem
        unpacked_path = temp_dir_path / extaction_dir_name
        print(f'  Tar detected; extracting {par_value} to {unpacked_path}', flush=True)

        with tarfile.open(par_value, 'r') as open_tar:
          members = open_tar.getmembers()
          root_dirs = [member
            for member in members 
            if member.isdir() and member.name != '.' and '/' not in member.name]
          # if there is only one root_dir (and there are files in that directory)
          # strip that directory name from the destination folder
          if len(root_dirs) == 1:
            for mem in members: 
              mem.path = Path(*Path(mem.path).parts[1:])
          members_to_move = [mem for mem in members if mem.path != Path('.')]
          open_tar.extractall(unpacked_path, members=members_to_move)
        new_values.append(unpacked_path)
      elif Path(par_value).is_file() and is_gz_file(par_value):
        # Remove extension (if it exists)
        extaction_file_name = Path(par_value.stem)
        unpacked_path = temp_dir_path / extaction_file_name
        print(f'  Gzip detected; extracting {par_value} to {unpacked_path}', flush=True)
        
        with gzip.open(par_value, 'rb') as f_in:
          with open(unpacked_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        new_values.append(unpacked_path)
      else:
        new_values.append(par_value)
    
    # unlist if need be
    if not is_multiple:
      new_values = new_values[0]

    # replace value
    par[par_name] = new_values
  # end ungzipping
  print("", flush=True)

  print("Grouping R1/R2 input files into pairs", flush=True)
  input_grouped = {}
  for file in par['readFilesIn']:
    key = re.search(fastqgz_regex, file.name).group(4)
    if key not in input_grouped:
      input_grouped[key] = []
    input_grouped[key].append(str(file))
  par['readFilesIn'] = [ ','.join(val) for val in input_grouped.values() ]
  print("", flush=True)

  print(">> Constructing command", flush=True)
  out_tmp_dir = temp_dir_path / "run"
  cmd_args = [
    # "/opt/cellranger-7.0.1/lib/bin/STAR",
    "STAR",
    "--runMode",
    "alignReads",
    "--outTmpDir",
    out_tmp_dir
  ]
  if 'cpus' in meta and meta['cpus']:
    cmd_args.extend(["--runThreadN", meta['cpus']])

  for name, value in par.items():
    if value is not None:
        if isinstance(value, list):
          cmd_args.extend(["--" + name] + value)
        else:
          cmd_args.extend(["--" + name, value])
  print("", flush=True)

  print(">> Running STAR with command:", flush=True)
  print("+ " + ' '.join([str(x) for x in cmd_args]), flush=True)
  print("", flush=True)

  subprocess.run(
    cmd_args,
    check=True
  )