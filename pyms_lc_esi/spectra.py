#!/usr/bin/env python3
#
#  spectra.py
"""
Mass spectrum handling.
"""
#
#  Copyright © 2020-2023 Dominic Davis-Foster <dominic@davis-foster.co.uk>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
#  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
#  OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
#  OR OTHER DEALINGS IN THE SOFTWARE.
#

# stdlib
# from typing import List, NamedTuple, Tuple

# 3rd party
# import numpy
# import pandas
from chemistry_tools.formulae import IsotopeDistribution
# from mathematical.data_frames import set_display_options
# from matplotlib.axes import Axes
# from matplotlib.container import BarContainer
# from pyms import Peak
# from pyms.IntensityMatrix import BaseIntensityMatrix
# from pyms.Spectrum import CompositeMassSpectrum, MassSpectrum, normalize_mass_spec
from pyms.Spectrum import MassSpectrum  # type: ignore[import]

__all__ = ["iso_dist_2_mass_spec"]

# set_display_options()


def iso_dist_2_mass_spec(
		iso_dist: IsotopeDistribution,
		min_abundance: float = 0,
		) -> MassSpectrum:
	"""
	Returns the Mass Spectrum representation of the given isotope distribution.

	:param iso_dist:
	:param min_abundance: Ignore isotopologues whose (absolute) abundance is below this threshold.
		By default isotopologues with zero abundance are excluded.
	:no-default min_abundance:
	"""

	iso_df = iso_dist.as_dataframe(format_percentage=False)
	iso_df = iso_df.astype({"Mass": float, "Abundance": float, "Relative Abundance": float})
	iso_df = iso_df[iso_df["Abundance"] > min_abundance]
	iso_df.sort_values(by=["Mass"], axis=0, ascending=True, inplace=True)
	iso_df.reset_index(inplace=True, drop=True)

	# print(iso_df)

	# return MassSpectrum(iso_df["Mass"], iso_df["Abundance"])
	return MassSpectrum(iso_df["Mass"], iso_df["Relative Abundance"])


# class LabelledMassSpectrum(NamedTuple):
# 	spectrum: MassSpectrum
# 	label: str

# def plot_isotope_distribution(
# 		ax: Axes,
# 		experimental: MassSpectrum,
# 		predicted: MassSpectrum,
# 		label: str = "predicted",
# 		mass_tol: float = 0.075
# 		) -> Tuple[BarContainer, BarContainer]:

# 	# TODO: for experimental, only plot peaks in common with predicted, even if more in the plot

# 	experimental_plot = plot_mass_spec(
# 			ax,
# 			normalize_mass_spec(experimental, max_intensity=100.0),
# 			label="experimental",
# 			width=0.025,
# 			color="C0",
# 			edgecolor="white",
# 			)

# 	predicted_plot = plot_mass_spec(
# 			ax,
# 			normalize_mass_spec(predicted, max_intensity=100.0),
# 			label=label,
# 			width=mass_tol,
# 			facecolor="None",
# 			edgecolor="C1",
# 			)

# 	return experimental_plot, predicted_plot

# #
# #
# # def plot_isotope_distribution(
# # 		ax: Axes,
# # 		experimental: MassSpectrum,
# # 		*predicted: LabelledMassSpectrum,
# # 		):
# #
# # 	min_mass = experimental.min_mass
# # 	max_mass = experimental.max_mass
# #
# # 	plot_mass_spec(
# # 			ax,
# # 			normalize_mass_spec(experimental, max_intensity=100.0),
# # 			label="experimental",
# # 			width=0.025,
# # 			color="C0",
# # 			edgecolor="white",
# # 			)
# #
# # 	for spec, label in predicted:
# # 		plot_mass_spec(
# # 				ax,
# # 				normalize_mass_spec(spec, max_intensity=100.0),
# # 				label=label,
# # 				width=0.1,
# # 				facecolor="None",
# # 				edgecolor="C1",
# # 				)
# #
# # 		min_mass = min((min_mass, spec.min_mass))
# # 		max_mass = max((max_mass, spec.max_mass))
# #
# # 	ax.set_xlim(min_mass-1, max_mass+1)

# def compare_iso_dist(predicted: MassSpectrum, experimental: MassSpectrum) -> pandas.DataFrame:
# 	predicted = normalize_mass_spec(predicted, max_intensity=100.0)
# 	experimental = normalize_mass_spec(experimental, max_intensity=100.0)

# 	df = pandas.DataFrame(
# 			data={
# 					"Predicted Masses": predicted.mass_list,
# 					"Predicted Intensities": predicted.mass_spec,
# 					}
# 			)
# 	df["Experimental Masses"] = numpy.nan
# 	df["Experimental Intensities"] = numpy.nan

# 	# Remove masses with zero normalized intensity
# 	df = df[df["Predicted Intensities"] > 0.0]
# 	df.sort_values(by=["Predicted Masses"], axis=0, ascending=True, inplace=True)
# 	df.reset_index(inplace=True, drop=True)

# 	extra_peaks = {}

# 	for mass, intensity in zip(experimental.mass_list, experimental.mass_spec):
# 		idx = (numpy.abs(df["Predicted Masses"] - mass)).argmin()
# 		if numpy.isclose(df["Predicted Masses"][idx], mass, atol=0.05):
# 			if numpy.isnan(df.at[idx, "Experimental Masses"]):
# 				# No existing value
# 				df.at[idx, "Experimental Masses"] = mass
# 				df.at[idx, "Experimental Intensities"] = intensity
# 			else:
# 				# Don't clobber existing value
# 				if df.at[idx, "Experimental Intensities"] > intensity:
# 					extra_peaks[mass] = intensity
# 				else:
# 					extra_peaks[df.at[idx, "Experimental Masses"]] = df.at[idx, "Experimental Intensities"]
# 					df.at[idx, "Experimental Masses"] = mass
# 					df.at[idx, "Experimental Intensities"] = intensity
# 		else:
# 			extra_peaks[mass] = intensity

# 	df.sort_values(by=["Predicted Masses"], axis=0, ascending=True, inplace=True)
# 	df.reset_index(inplace=True, drop=True)

# 	df["Mass Δ"] = df.apply(func=df_delta, args=["Predicted Masses", "Experimental Masses"], axis=1)
# 	df["Rel. Mass Δ"] = df.apply(func=df_delta_relative, args=["Predicted Masses", "Experimental Masses"], axis=1)
# 	df["Intensity Δ"] = df.apply(func=df_delta, args=["Predicted Intensities", "Experimental Intensities"], axis=1)
# 	df["Rel. Intensity Δ"] = df.apply(
# 			func=df_delta_relative, args=["Predicted Intensities", "Experimental Intensities"], axis=1
# 			)

# 	df.sort_values(by=["Predicted Intensities"], axis=0, ascending=False, inplace=True)
# 	df.reset_index(inplace=True, drop=True)

# 	# print(df)
# 	# print("Average Mass Difference:", numpy.nanmean(df["Mass Δ"]))
# 	# print("Average Intensity Difference:", numpy.nanmean(df["Intensity Δ"]))
# 	# pprint(extra_peaks)

# 	return df

# def df_delta(row: pandas.Series, left_column: str, right_column: str) -> float:
# 	"""
# 	Calculate the difference between values in the two columns for each row of a
# 	:class:`data frame <pandas.DataFrame>`.

# 	Do not call this function directly; use it with
# 	:meth:`df.apply() <pandas.DataFrame.apply>` instead:

# 	.. code-block:: python

# 		data_frame["Mean"] = data_frame.apply(
# 				func=df_delta,
# 				args=["Bob", "Alice"],
# 				axis=1,
# 				)

# 	:param row: Row of the data frame.
# 	:param left_column:
# 	:param right_column:

# 	:return: The difference between ``left_column`` and ``right_column``.
# 	"""  # noqa D400

# 	return row[left_column] - row[right_column]

# def df_delta_relative(row: pandas.Series, left_column: str, right_column: str) -> float:
# 	"""
# 	Calculate the relative difference between values in the two columns for each row of a
# 	:class:`data frame <pandas.DataFrame>`.

# 	Do not call this function directly; use it with
# 	:meth:`df.apply() <pandas.DataFrame.apply>` instead:

# 	.. code-block:: python

# 		data_frame["Mean"] = data_frame.apply(
# 				func=df_delta_relative,
# 				args=["Bob", "Alice"],
# 				axis=1,
# 				)

# 	:param row: Row of the data frame.
# 	:param left_column:
# 	:param right_column:

# 	:return: The relative difference between ``left_column`` and ``right_column``.
# 	"""  # noqa D400

# 	right = row[right_column]
# 	if right:
# 		return (row[left_column] - right) / right
# 	else:
# 		return float("inf")

# def find_best_peak(peak_list: List[Peak.Peak],
# 					reference: MassSpectrum) -> List[Tuple[Peak.Peak, pandas.DataFrame]]:
# 	results = [(peak, compare_iso_dist(reference, peak.mass_spectrum)) for peak in peak_list]

# 	def sort_func(item: Tuple[Peak.Peak, pandas.DataFrame]):
# 		df = item[1]
# 		return 100.0 / df.at[0, "Rel. Mass Δ"
# 								]  # ,  100 / df.at[1, "Rel. Intensity Δ"]  # , 100.0 / df.at[1, "Rel. Mass Δ"]

# 	results.sort(key=sort_func)

# 	# TODO: break ties. However, sorting a list of DataFrames is not stable,
# 	# so probably need a custom sort function that uses
# 	# https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.equals.html

# 	return results

# def composite_spectrum_from_peak(peak: Peak.Peak, im: BaseIntensityMatrix) -> "CompositeMassSpectrum":
# 	bounds = peak.bounds
# 	if not bounds or not bounds[0] or not bounds[1]:
# 		raise ValueError("The peak has no bounds set!")

# 	left_bound, apex, right_bound = peak.bounds
# 	left_bound = apex - left_bound
# 	right_bound = apex + right_bound

# 	print(peak.bounds)
# 	spectra: List[MassSpectrum] = []
# 	for idx in range(left_bound, right_bound + 1):
# 		print(idx, im.get_ms_at_index(idx))
# 		spectra.append(im.get_ms_at_index(idx))

# 	# print(reduce(operator.add, spectra))
# 	return CompositeMassSpectrum.from_spectra(spectra)
