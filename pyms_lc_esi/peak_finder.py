# stdlib
from itertools import chain
from typing import Iterable, Iterator, Tuple

# 3rd party
from chemistry_tools.formulae import Formula
from domdf_python_tools.words import word_join
from pyms.BillerBiemann import get_maxima_indices, num_ions_threshold
from pyms.eic import ExtractedIntensityMatrix, build_extracted_intensity_matrix
from pyms.IntensityMatrix import IntensityMatrix
from pyms.Noise.Analysis import window_analyzer
from pyms.Peak import Peak

# this package
from pyms_lc_esi.adducts import Adduct, get_adduct_spectra

__all__ = ["make_im_for_adducts", "peak_finder", "peaks_from_maxima", "sum_area"]


def sum_area(
		apex_index: int,
		e_im: ExtractedIntensityMatrix,
		) -> Tuple[float, int, int]:
	"""
	Returns the area and absolute bounds (as scans) for the peak with the apex at ``apex_index``.

	.. TODO:: explain how bounds are determined

	:param apex_index: The scan index of the apex of the peak.
	:param e_im:
	"""

	# print(peak, peak.rt / 60)

	# print(f"Apexes at {apex_index}")
	# print("Intensity at apex:")

	apex_intensity = sum(e_im.intensity_array[apex_index])
	rhs_area = lhs_area = last_intensity = apex_intensity
	left_bound = right_bound = apex_index
	bound_area_tolerance = 0.0005 / 2  # half of 0.05 %

	for idx_offset, scan in enumerate(e_im.intensity_array[apex_index + 1:]):
		scan_intensity = sum(scan)

		if scan_intensity <= last_intensity and scan_intensity > (rhs_area * bound_area_tolerance):
			last_intensity = scan_intensity
			rhs_area += scan_intensity
			right_bound = idx_offset + apex_index + 1
		else:
			break

	last_intensity = apex_intensity

	for idx_offset, scan in enumerate(reversed(e_im.intensity_array[:apex_index])):
		scan_intensity = sum(scan)

		if scan_intensity <= last_intensity and scan_intensity > (lhs_area * bound_area_tolerance):
			last_intensity = scan_intensity
			lhs_area += scan_intensity
			left_bound = (apex_index - idx_offset) - 1
		else:
			break

	area = lhs_area + rhs_area - apex_intensity  # apex intensity was counted for each half

	# print("Peak area:")
	# print(area)

	# print(f"right bound: {right_bound}", eic.get_time_at_index(right_bound) / 60)
	# print(f"left bound: {left_bound}", eic.get_time_at_index(left_bound) / 60)
	return area, left_bound, right_bound


def peaks_from_maxima(e_im: ExtractedIntensityMatrix, points=3):
	"""

	:param e_im:
	:param points:
	"""

	# TODO: combine close peaks, via scans param.
	intensity_list = []
	peak_list = []

	for row in e_im._intensity_array:
		intensity_list.append(sum(row))

	for apex_idx in get_maxima_indices(intensity_list, points=points):
		rt = e_im.get_time_at_index(apex_idx)
		ms = e_im.get_ms_at_index(apex_idx)
		peak = Peak(rt, ms)
		peak.bounds = (0, apex_idx, 0)

		peak_list.append(peak)

	return peak_list


def make_im_for_adducts(
		im: IntensityMatrix,
		analyte: Formula,
		adducts: Iterable[Adduct],
		left_bound: float = 0.1,
		right_bound: float = 0.1,
		):
	"""
	Conxtructs a :class:`pyms.eic.ExtractedIntensityMatrix` for the given adducts of the analyte.

	:param im:
	:param analyte:
	:param adducts:
	:param left_bound:
	:param right_bound:
	"""

	# Compile a list of masses for the adducts
	spectra = get_adduct_spectra(analyte, adducts)

	print(spectra)
	all_masses = sorted(set(chain.from_iterable(spectrum.mass_list for spectrum in spectra.values())))

	print(f"Constructing ExtractedIntensityMatrix for m/z {word_join(map(str, all_masses))}")

	# Construct the extracted intensity matrix for the adducts
	return build_extracted_intensity_matrix(im, masses=all_masses, left_bound=left_bound, right_bound=right_bound)


def peak_finder(e_im: ExtractedIntensityMatrix, points: int = 3) -> Iterator[Peak]:

	# Find peaks
	# peaks = BillerBiemann(e_im, points=8, scans=3)
	peaks = peaks_from_maxima(e_im, points=points)

	eic = e_im.eic

	# Filter small peaks from peak list
	noise_level = window_analyzer(eic)
	print(f"Filtering peaks with fewer than {2}/{len(e_im.mass_list)} masses.")
	print("noise_level:", noise_level)
	filtered_peak_list = num_ions_threshold(peaks, n=2, cutoff=noise_level)
	# filtered_peak_list = num_ions_threshold(peaks, n=1, cutoff=noise_level)
	# filtered_peak_list = num_ions_threshold(peaks, n=3, cutoff=5000)
	# filtered_peak_list = peaks

	for peak in filtered_peak_list[::-1]:
		apex_index = eic.get_index_at_time(peak.rt)

		# Estimate peak area
		peak.area, left_bound, right_bound = sum_area(apex_index, e_im)

		# Assign bounds to peak as offsets.
		peak.bounds = [apex_index - left_bound, apex_index, right_bound - apex_index]

		if (right_bound - left_bound) > 3 and peak.area > 1000:
			# TODO: make 1000 dependent on data

			yield peak


# def fill_peak(eic: ExtractedIonChromatogram, peak: Peak, ax: Optional[Axes] = None) -> PolyCollection:
# 	if ax is None:
# 		# 3rd party
# 		import matplotlib.pyplot
# 		ax = matplotlib.pyplot.gca()

# 	left_bound, apex, right_bound = peak.bounds
# 	left_bound = apex - left_bound
# 	right_bound = apex + right_bound

# 	# Fill in the peak
# 	time_range = [eic.get_time_at_index(idx) / 60 for idx in (range(left_bound, right_bound + 1))]
# 	intensities = eic.intensity_array.tolist()[left_bound:right_bound + 1]
# 	return ax.fill_between(time_range, intensities)
