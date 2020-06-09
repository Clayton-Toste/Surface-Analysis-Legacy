from surface_analysis import Surface_Analysis
from graph import Graph
from sys import argv
from pathlib import Path
from os import path, mkdir, chdir, remove
from shutil import copyfile
from time import time, ctime, strftime
from datetime import timedelta
from multiprocessing import Pool, cpu_count
import logging, subprocess, typing


def analyze(quiet: bool=False, reset: bool=False, multiprocess: bool=False, pymol_filter: bool= False, protein_codes: set={}, protein_files: set={}, graphing: set={}, logger: logging.Logger=None) -> None:
	
	for arg in graphing:
		if arg not in {'distances_outer', 'distances_inner', 'touching_outer', 'touching_inner', 'sizes_outer', 'sizes_inner', 'save', 'show', 'stats', 'use_radius', 'log_scale', 'remove_outliers'}:
			raise ValueError("Graphing argument, '"+arg+"', is invalid")

	proteins = set()

	for protein in protein_codes:
		proteins.add((protein.lower(), None))
	for file in protein_files:
		proteins.add((file.stem.lower(), file))
		if not file.exists():
			raise ValueError("File, '"+str(file)+"', does not exist")


	start = time()
	if not proteins:
		raise ValueError('No protein argument.')

	try:
		chdir('proteins')
	except:
		try:
			mkdir('proteins')
			chdir('proteins')
		except:
			raise OSError("Could not access 'proteins' folder.")
	if not logger:
		logger = logging.getLogger('')
		logger.setLevel(logging.DEBUG)
		formatter = logging.Formatter('%(asctime)s -- %(offset)s -- %(code)s: %(message)s', datefmt='%H:%M:%S')
		handlers = [logging.FileHandler('../logs/{0}.log'.format(strftime('%m.%d.%Y.%H.%M.%S')))]
		if not quiet:
			handlers.append(logging.StreamHandler())
		for handler in handlers:
			handler.setFormatter(formatter)
		logger.handlers = handlers

	end = time()
	logger.debug('Initialization complete, beginning processing of proteins with codes: '+repr(list(protein[0] for protein in proteins))[1:-1]+((' with '+('graphing arguments: '+repr(graphing)[1:-1] if graphing else 'no graphing arguments')) if graphing else "") +'.', extra={'offset': timedelta(seconds=end-start), 'code': 'Main'})

	for protein, file in proteins:
		try:
			chdir(protein)
		except:
			try:
				mkdir(protein)
				chdir(protein)
			except:
				raise OSError("Could not access '"+protein+"' folder.")
		if reset:
			try:
				remove('.finished')
			except:
				pass
		elif path.isfile(protein+'.wrl'):
			end = time()
			logger.debug('Skipping fetching of model for protein.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
			chdir('..')
			break
		if not path.isfile('.file'):
			try:
				end = time()
				logger.debug('Building model.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
				if file:
					copyfile(file, file.name)
					open(".file", 'w').close()
				subprocess.call(['../../PyMOL/python.exe', '../../fetch.py', file.name if path.isfile('.file') else protein, str(not path.isfile('.file'))], stdout=subprocess.PIPE)
				end = time()
				logger.debug('Successfully built model.', extra={'offset': timedelta(seconds=end-start), 'code': protein})
			except:
				raise SystemError('Could not run fetch script.')
		chdir('..')

	end = time()
	logger.debug('Beginning surface analysis.', extra={'offset': timedelta(seconds=end-start), 'code': 'Main'})

	def process_protein(protein: str) -> None:
		if not path.isfile(protein[0]+'/.finished'):
			Surface_Analysis(protein[0], logger, start)
			return
		end = time()
		logger.debug('Skipping surface analysis for protein because it has already been processed.', extra={'offset': timedelta(seconds=end-start), 'code': protein[0]})

	if multiprocess:
		with Pool(cpu_count() if cpu_count() <= len(proteins) else len(proteins)) as p:
			p.map(process_protein, proteins)
	else:
		for protein in proteins:
			process_protein(protein)
		
	end = time()
	logger.debug('Surface analysis complete.'+' Beginning graphing.' if graphing else '', extra={'offset': timedelta(seconds=end-start), 'code': 'Main'})

	if graphing:
		Graph(set(protein for protein, file in proteins), graphing, logger, start)


if __name__ == '__main__':
	quiet = False
	reset = False
	multiprocess = False
	pymol_filter = False
	protein_codes = set()
	protein_files = set()
	graphing = set()
	mode = 0
	for arg in argv[1:]:
		if arg == '-h':
			print('''
				Usage: Main.py [arg] ...

				-h     : print this help message and exit
				-p     : add a protein to analyze by it's PDB id or the file name of an already processed protein, can have mutiple arguments
				-f     : add a protein to analyze from a file path, can have mutiple arguments
				-g     : add graphing options, can have mutiple arguments
				-r     : force data to recalulated even if already existing
				-q     : run quietly with no output to the console
				-i     : add a filter using python script in \\filters folder by filename (no extension) 

				*Graphing options:
				
				save                : save graphs to \\images folder
				show                : display graphs in an window
				remove_outliers     : remove outliers from data when graphing
				stats               : save statistics for graph to an fgile in \\statistics folder
				use_radius          : use radius of appoxmate circles to estimate area
				log_scale           : use log scale when ploting
				[type]_[location]   : the data to be graphed
				
				Types:
				distances  : distances between patchs the with the classification
				touching   : distances between patchs that are touching
				sizes      : sizes of patches
				
				Locations:
				Outer      : outer surface of protein
				Inner      : inner surface of protein

				*When graphing on windows Xming must be installed and running
			''')
			exit()
		if arg == '-q':
			quiet = True
			mode = 0
			continue
		if arg == '-r':
			reset = True
			mode = 0
			continue
		if arg == '-m':
			multiprocess = True
			mode = 0
			continue
		if arg == '-i':
			pymol_filter = True
			mode = 0
			continue
		if arg == '-p':
			mode = 1
			continue
		if arg == '-f':
			mode = 2
			continue
		if arg == '-g':
			mode = 3
			continue
		if mode == 1:
			protein_codes.add(arg)
			continue
		if mode == 2:
			file = arg
			if ':' == arg[1].lower():
				file = Path('/mnt/'+arg[1]+arg[2:])
			else:
				file = Path(arg).resolve()
			protein_files.add(file)
			continue
		if mode == 3:
			graphing.add(arg)
			continue
		raise ValueError('Argument "'+arg+'" is invalid.')
	analyze(quiet=quiet, reset=reset, multiprocess=multiprocess, protein_codes=protein_codes, protein_files=protein_files, graphing=graphing, pymol_filter=pymol_filter)
