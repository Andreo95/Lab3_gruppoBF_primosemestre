import numpy as np

def maketab( *columns , errors=None, precision=3):

	vals = np.asarray(columns).T
	cols = len(vals.T)
	precision = precision * np.ones(cols)
	if errors in ('all', 'any'): errors = range(1, cols, 2)
	elif errors in ('yes','y', 1): errors = [1]
	
	tab = (
		'\\begin{table}\n'+
		'\t\\begin{tabular}{' +
		'*{{{0}}}'.format(cols if errors is None else cols-len(errors) ) +
		'{S}} \n\n' +
		'\t\t\\midrule \n'
	)

	for row in vals:
		rows = enumerate(row, start=1)
		tab += '\t'
		for pos,v in rows:
			num = v if np.isfinite(v) else '-'
			space = '&' if pos < cols else '\\\\ \n'
			err = ''
			prec = None
			
			if errors is not None and pos in errors:
				prec = np.floor(np.log10(row[pos]))
				if row[pos]/10**prec < 2.5:
					prec -= 1
				err = '({0:.0f})'.format(round(row[pos]/10**prec))
				if next(rows, (-1, None))[0] >= cols : space = '\\\\ \n'
				num = round(num, int(-prec))
				
			tab += '\t{0:.{digits:.0f}f} {1}\t{2}'.format(num, err, space, digits = precision[pos-1] if prec is None else max(0, -prec) )
		
	tab += '\t\\end{tabular} \n' + '\t\\caption{some caption} \n' + '\t\\label{t:somelabel} \n' + '\\end{table}'
	
	return tab
