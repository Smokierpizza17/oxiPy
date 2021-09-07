

# oxiPy
Python program for computing oxidisation states for (anorganic) chemical molecules

can be used from CLI:
```
> py oxi.py "Na(NO3)"
1   5 -2  
Na (N  O₃)
```
or as an import in python:
```python
>>> import oxi
>>> print(oxi.wrapper("Na(NO3)"))
1   5 -2  
Na (N  O₃)
```
## input
Chemical formulas are written without spaces and with plain numbers instead of subscript. Overall charge is written after a space. Subgroups (polyatomic ions that can be solved from a lookup table) should be written in brackets.

### examples
|molecule |formula |input string |
--- | --- | ---
|sodium ion |Na⁺ | ```"Na +"```, ```"Na 1+"``` |
|copper(II) ion |Cu²⁺ | ```Cu 2+``` |
|carbonate |CO₃²⁻ | ```"CO3 2-"``` |
|sodium carbonate |Na₂CO₃ |```"Na2(CO3)"``` |
|copper(II) nitrate | Cu(NO₃)₂ |```Cu(NO3)2``` |

## structure
There are three important functions: ```getOxiNumbers```, ```interpretElements```, and ```printResult```.

```interpretElements``` splits up the provided string and returns the formula as ```([[elementSymbol, count, None], [subgroupName, count, None], ...], overallCharge)```

```python
interpretElements("Na(NO3)")  # --> ([['Na', 1, None], ('NO3', 1, None)], 0)
```

```getOxiNumbers``` does everything that ```interpretElements``` does and computes all the oxidisation states. For subgroups, they are shown as sets with their constituent atoms and respective oxidisation states.

```python
getOxiNumbers("Na(NO3)") # --> ([['Na', 1, 1], [(['N', 1, 5], ['O', 3, -2]), 1, -1]], 0)
```

```printResult``` takes an array like the one returned by ```getOxiNumbers``` and pretty-prints the result with superscript charge and subscript indexes, as well as the oxidisation states aligned above the elements.
```python
# tuple returned by getOxiNumbers has to be provided as separate arguments with *
printResult(*getOxiNumbers("Na(NO3)")) # --> 1   5 -2
                                             Na (N  O₃)
```

## found an error?
Please create an issue with the label ```algorithmError```!
