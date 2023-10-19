#!/usr/bin/env python3
#
#  adducts.py
"""
Represents LC-MS adducts.
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
from typing import Dict, Iterable, Literal, Union

# 3rd party
import attr
from chemistry_tools.elements import ELEMENTS
from chemistry_tools.formulae import Formula
from pyms.Spectrum import MassSpectrum

# this package
from pyms_lc_esi.spectra import iso_dist_2_mass_spec

__all__ = [
		"Adduct",
		"plus_h",
		"plus_sodium",
		"get_adduct_spectra",
		]


def _formula_converter(formula: Union[Dict[str, int], Formula, str]) -> Formula:
	if isinstance(formula, Formula):
		return formula
	elif isinstance(formula, Dict):
		return Formula.from_kwargs(**formula)
	elif isinstance(formula, str):
		if formula not in ELEMENTS:
			raise ValueError(f"{formula!r} is not an Element")
		return Formula.from_kwargs(**{formula: 1})
	else:
		raise TypeError(f"Unsupported type for formula: {type(formula)}")


def _operation_validator(instance, attribute, value: Literal["add", "sub"]) -> Literal["add", "sub"]:
	if not isinstance(value, str):
		raise TypeError(f"operation must be a string, not {type(value)}")

	value = value.lower().strip()
	if value not in {"add", "sub"}:
		raise ValueError(f"Unsupported operation {value}")

	return value


@attr.s
class Adduct:
	"""
	Create adducts of formulae.
	"""

	#: The name of the adduct, e.g. ``'[%s + H]⁺'``.
	name: str = attr.ib(converter=str)

	formula: Union[Dict[str, int], Formula, str] = attr.ib(converter=_formula_converter)
	"""
	Either:

	* a :class:`chemistry_tools.formulae.formula.Formula` representing the elements and their quantities to
	  be added/removed to create the adduct;
	* a dict representing the same; or
	* a string giving the symbol of an element, of which a single atom is added/removed to create the adduct.
	"""

	#:
	operation: Literal["add", "sub"] = attr.ib(validator=_operation_validator, default="add")

	def __call__(self, formula: Formula) -> Formula:
		"""
		Create the adduct.

		:param formula: The base formula.

		:returns: The formula of the adduct.
		"""

		if self.operation == "add":
			return formula + self.formula
		elif self.operation == "sub":
			return formula - self.formula

	def __format__(self, format_spec):
		return self.name % format_spec

	def __mod__(self, other):
		return self.name % other


# Some common adducts
plus_h = Adduct("[%s + H]⁺", 'H')
plus_sodium = Adduct("[%s + Na]⁺", "Na")


def get_adduct_spectra(
		formula: Formula,
		adducts: Iterable[Adduct],
		) -> Dict[str, MassSpectrum]:
	"""
	Returns a dictionary mapping adducts to mass spectra, for the given adducts of ``formula``.

	:param formula:
	:param adducts:
	"""

	# TODO: Add parameter to this and all calling functions to control cutoff intensity.

	spectra = {}
	print(formula)

	for adduct in adducts:
		formula_as_adduct = adduct(formula)
		# print(f"{adduct % 'M'}:", formula_as_adduct)
		# print("Average Mass:", formula_as_adduct.mass)
		# print("Exact Mass:", formula_as_adduct.exact_mass)

		# print()

		# print(f"Isotope distribution of {adduct % 'Diphenylamine'}:")
		# print(formula_as_adduct.isotope_distribution())
		dpa_h_mass_spec = iso_dist_2_mass_spec(formula_as_adduct.isotope_distribution(), 0.001)

		spectra[adduct % 'M'] = dpa_h_mass_spec

	return spectra
