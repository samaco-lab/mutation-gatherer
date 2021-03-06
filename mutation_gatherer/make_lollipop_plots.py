#derived from Parashar Dhapola's https://gist.github.com/parashardhapola/e095e3ebd068c90e7089b53e46b8e0bc

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib
from itertools import repeat

class GenomicLollipopPlot(object):
	def __init__(self, exon_intervals, marker_pos=[], marker_heights=[], marker_colors=[],
				 marker_size=100, marker_weight=1.5, exon_color="black", intron_color="grey",
				 intron_weight=2, intron_style='-', bar_color='cornflowerblue', bg_color="white"):
		self.exonIntervals = exon_intervals
		self.markerPositions = marker_pos
		self.markerHeights = marker_heights
		self.markerColors = marker_colors
		self.markerSize = marker_size
		self.MarkerWeight = marker_weight
		self.exonColor = exon_color
		self.intronColor = intron_color
		self.intronWeight = intron_weight
		self.intronStyle = intron_style
		self.barColor= bar_color
		self.bgColor = bg_color
		self.markerDefaultColor = 'grey'
		self.numExons = len(self.exonIntervals)
		self.totalSpan = self.exonIntervals[-1][1] - self.exonIntervals[0][0]
		self.minExonLen = self.totalSpan*0.005
		self.ylims = {'exon_max': 2, 'exon_min':1}
		self.figure, self.canvas = plt.subplots(figsize=(15,1.5))
		#self.canvas.set_axis_bgcolor(self.bgColor)
		self._draw()

	def _set_limits(self):
		self.ylims['intron_max'] = self.ylims['exon_max']*0.9
		self.ylims['intron_min'] = (self.ylims['exon_max'] + self.ylims['exon_min'])/2.0
		self.ylims['bar_min'] = self.ylims['exon_max']+0.2
		self.ylims['bar_max'] = self.ylims['bar_min']+(self.ylims['exon_max']-self.ylims['exon_min'])/5.0
		
	
	def _transform_spans(self):
		span_lens = [x[1]-x[0] for x in self.exonIntervals]
		max_len = float(max(span_lens))
		transformed_intervals = []
		if max_len < self.minExonLen:
			span_ratios = [x/max_len for x in span_lens]
			expansion_factor = self.totalSpan*1e-11
			for i in range(1,10):
				ef = (2**i)*expansion_factor
				if max_len+ef > self.minExonLen:
					expansion_factor = ef
					break
			for i,j in zip(self.exonIntervals, span_ratios):
				mid = (i[0] + i[1])/2
				f = (expansion_factor*j)/2
				if mid+f - mid-f > self.minExonLen:
					transformed_intervals.append([mid-f, mid+f])
				else:
					transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
		else:
			for i in range(self.numExons):
				if span_lens[i] < self.minExonLen:
					mid = (self.exonIntervals[i][0] + self.exonIntervals[i][0])/2 
					transformed_intervals.append([mid-(self.minExonLen/2), mid+(self.minExonLen/2)])
				else:
					transformed_intervals.append(self.exonIntervals[i])
		self.exonIntervals = transformed_intervals[:]
		
	def _draw_exon(self, span):
		self.canvas.fill_between(span, self.ylims['exon_min'], self.ylims['exon_max'],
								 edgecolor=self.bgColor, facecolor=self.exonColor)
		return True
		
	def _draw_intron(self, span):
		mid = (span[0]+span[1])/2.0
		self.canvas.plot([span[0], mid], [self.ylims['intron_min'], self.ylims['intron_max']],
						 c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		self.canvas.plot([mid, span[1]], [self.ylims['intron_max'], self.ylims['intron_min']],
						 c=self.intronColor, lw=self.intronWeight, ls=self.intronStyle)
		return True
	
	def _draw_markers(self):
		if self.markerHeights == []:
			self.markerHeights = [self.ylims['exon_max']-self.ylims['exon_min'] for x in self.markerPositions]
		if self.markerColors == []:
			self.markerColors = [self.markerDefaultColor for x in self.markerPositions]           
		for p,h,c in zip(self.markerPositions, self.markerHeights, self.markerColors):
			self.canvas.plot((p, p), (self.ylims['bar_max'], self.ylims['bar_max']+h),
							 linestyle='-', color='black', linewidth=self.MarkerWeight, alpha=0.7)
			self.canvas.scatter(p, self.ylims['bar_max']+h+0.25, s=self.markerSize, marker='o', c=c,
								edgecolor=c, alpha=1)
		
	
	def _clean_axes(self):
		#self.canvas.set_ylim((self.ylims['exon_min'], self.ylims['bar_max']))
		#self.canvas.set_yticks([], [])
		self.canvas.get_xaxis().tick_top()
		self.canvas.tick_params(axis='x', direction='out')
		self.canvas.set_xticks([])
		for o in ["top", "bottom", "left", "right"]:
			self.canvas.spines[o].set_visible(False)
		min_pos = int(self.exonIntervals[0][0] - self.totalSpan * 0.1)
		if min_pos < 0:
			min_pos = 0
		max_pos = int(self.exonIntervals[-1][1] + self.totalSpan * 0.1)
		minortick_pos = [x for x in range(min_pos, max_pos, int((max_pos-min_pos)/20))][1:]
		for i in minortick_pos:
			self.canvas.axvline(i, alpha=0.2, c='black', ls='--')
		self.canvas.text(minortick_pos[0], self.ylims['exon_min']-0.5,
						 minortick_pos[0], fontsize=8, ha='center')
		self.canvas.text(minortick_pos[-1], self.ylims['exon_min']-0.5,
						 minortick_pos[-1], fontsize=8, ha='center')
		self.canvas.set_xlim(minortick_pos[0]-(minortick_pos[1]-minortick_pos[0]),
							 minortick_pos[-1]+(minortick_pos[-1]-minortick_pos[-2]))
		
	def _draw(self):
		self._set_limits()
		self._transform_spans()
		for i in range(self.numExons):
			if i > 0:
				self._draw_intron([self.exonIntervals[i-1][1], self.exonIntervals[i][0]])
			self._draw_exon(self.exonIntervals[i])
		self.canvas.fill_between([self.exonIntervals[0][0], self.exonIntervals[-1][1]],
								  self.ylims['bar_min'], self.ylims['bar_max'],
								  edgecolor=self.bgColor, facecolor=self.barColor)
		self._draw_markers()
		self._clean_axes()
	
	def show(self):
		plt.show()

def test():
	#exons positions
	exon_pos = [97543299,97544702], [97547885,97548026], [97564044,97564188], [97658624,97658804], [97700407,97700550], [97770814,97770934], [97771732,97771853], [97839116,97839200], [97847948,97848017], [97915614,97915779], [97981281,97981497], [98015115,98015300], [98039315,98039526], [98058773,98058943], [98060614,98060722], [98144650,98144738], [98157272,98157354], [98164906,98165103], [98187065,98187227], [98205947,98206035], [98293669,98293752], [98348819,98348930], [98386439,98386615]
	#marker positions
	marker_pos = [97947885, 98247485]
	gene = GenomicLollipopPlot(exon_pos, marker_pos)
	gene.show()    	


def obtain_markers():
	pass

def cdna_lollipop_plot(exon_coordinates, variant_coordinates=[10,100,1000,10000,15000]):
	fig = plt.figure(figsize=(50,1.5))
	ax = fig.add_subplot()
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['left'].set_visible(False)	

	x1 = 0
	for start, end in exon_coordinates:
		size = end - start
		x2 = x1 + size
		exon = matplotlib.patches.Rectangle((x1,0),x2,50, edgecolor = 'red', facecolor = 'paleturquoise')
		ax.add_patch(exon)
		x1 = x2

	xlim1 = -50
	xlim2 = x2 + 100

	plt.xlim([xlim1,xlim2])
	plt.ylim([-10,110])

	stems = []
	stems.extend(repeat(200,len(variant_coordinates)))
	plt.stem(variant_coordinates, stems, bottom = 50, linefmt = 'grey', use_line_collection = False, markerfmt = 'D')

	plt.show()


def main():
	pass
		
if  __name__ == "__main__":
	main()