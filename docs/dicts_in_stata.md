# Dicts in Stata

Dictionaries, also known as associative arrays, maps, or symbol tables have an unusual syntax in Stata, compared to e.g. Python.

## Python Equivalency

```python
A = dict()
```
```stata
A = asarray_create([keytype,keydim,minsize,minratio,maxratio])

```