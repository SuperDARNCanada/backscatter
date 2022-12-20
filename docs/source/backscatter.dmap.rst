backscatter.dmap package
========================
Python based library for reading and writing SuperDARN data in the dmap format

This project attempts to recreate the DMAP file encoding and decoding from RST. This version requires no dependencies on RST and can easily be imported to your Python scripts. The only dependency this project has is Numpy, as the package is both compatible with Numpy arrays and standard Python lists. This project is much more robust in error handling than the C version so it can be used to parse radar data from sockets as well.

DMAP records have a packet style structure and stored as a binary string. Because each DMAP record stores its size, records can be concatenated together in a file or buffer and still be decoded. DMAP records have the following structure and data types:

====  ====  =========  ========  =============  =============
code  size  # scalers  # arrays  scalers        arrays 
====  ====  =========  ========  =============  =============
int   int   int        int       variable size  variable size
====  ====  =========  ========  =============  =============

Scalers and arrays also have a packet style structure. Scalers look like this:

=============  ====  ====
name           type  data
=============  ====  ====
variable size  char  type
=============  ====  ====

To parse a scaler, bytes are read as a string until a null terminator is encountered. This is name of the scaler. Then a single char is read which corresponds to the data type. A number of bytes corresponding to that type is then parsed as data.

Arrays look like this:

=============  ====  ============  ===============  =====================
name           type  # dimensions  dimensions       data
=============  ====  ============  ===============  =====================
variable size  char  int           variable # ints  total elements * type
=============  ====  ============  ===============  =====================

To parse an array, bytes are read as a string until a null terminator is encountered. This is name of the array. Then a single char is read which corresponds to the data type. An int is read to determine the number of array dimensions. For the number of array dimensions, ints are read to know what each array dimension is. Once all array dimensions are known, the total number of elements can be determined, then multiplied by size of parsed type, and the total determined bytes can be read into a contiguous piece of memory. The array can then be reshaped to the correct dimensions.

The following parsed type values will correspond to these data types:

====  ====  =====  ===  =====  ======  ======  ====  =====  ======  ====  =====   
0     1     2      3    4      8       9       10    16     17      18    19
====  ====  =====  ===  =====  ======  ======  ====  =====  ======  ====  =====
DMAP  CHAR  SHORT  INT  FLOAT  DOUBLE  STRING  LONG  UCHAR  USHORT  UINT  ULONG 
====  ====  =====  ===  =====  ======  ======  ====  =====  ======  ====  =====


Reading and decoding data
-------------------------
To work with the package there are simple to use functions that abstract the DMAP complexity for the user. To read and decode the data from any DMAP formatted file, use:
:code:`parse_dmap_format_from_file(filepath, raw_dmap=False)`

or if working on data from a buffer(such as data from a network socket) then use:
:code:`parse_dmap_format_from_stream(stream, raw_dmap=False)`


Both of these functions return a list of dictionaries that hold the data for each parsed DMAP record, or if :code:`raw_dmap` is set :code:`True` then it returns the underlying Raw_Dmap object that holds data while it is being parsed. Note that arrays are parsed using Numpy for speed and are of :code:`numpy.ndarray` type.

Encoding and writing data
-------------------------
The underlying algorithm for writing and encoding data to a DMAP file is slightly more complicated. When determining data types to be encoded, dmap writing will default to use the Python type, unless a user defined dictionary is supplied that overrides what type of data something is. This is neccessary to order to make files compatible with the C reading routines. For example in a .fitacf file, `gflg` is an array of 0s and 1s which Python will generate or store as ints. However, this needs to be written out as chars to stay compatible for the C version. The writing routines can also handle both Python lists and Numpy arrays. There is a function that abstracts this complexity for .iqdat,.rawacf,.fitacf, and .map types:
:code:`dicts_to_file(data_dicts, file_path, file_type='')`

This function takes a list of dictionaries to write, a file path to write to, and a file type to type override Python's default. Options are 'iqdat', 'rawacf', 'fitacf', and 'map'. 

To create a custom writing routine, a user can create a :code:`RawDmapWrite(data_dicts,file_path,ud_types={})` object and pass it your own user defined dictionary for types similar to above wrapper function. For type formats, refer to the :code:`struct` module, or if using Numpy arrays then refer to the Numpy :code:`dtypes`. When creating a user defined dictionary, key-value pairs should be in the form:
:code:`{
'name' : 'fmt', 
}`

where :code:`name` is the name of dictionary element to override, and :code:`fmt` is a type format listed in the :code:`struct` or :code:`dtypes` documentation.

Exceptions
----------
pydmap will raise :code:`DmapDataError` on failure. There are many places where the package can raise this exception so it is recommended wrapping all functions in a :code:`try` clause and handling the exception as you see fit.

Testing
-------
There are a set of basic unit tests if you make changes. These tests affirm that the package can open and parse common file types. It tests writing capabilites by encoding some data, writing it out, and then trying to parse it again. There is also a randomized test to test error handling of the parsing routines. This test takes a valid DMAP record and randomly corrupts 5% of the data to see if the parser throws any unhandled exceptions.

Unimplemented or untested features
----------------------------------
I have not tested recursive DMAP records. I'm not aware of any files that use DMAP records within records.

Differing byte endianess. I'm not aware of any files that were written using differing byte endianess. If this is an issue, it can be added

I have not full tested the backwards compatibility of all pydmap written files using the C data analysis programs. The few files I tested did work, but their may be some edge cases

Submodules
----------

.. toctree::

   backscatter.dmap.dmap

Module contents
---------------

.. automodule:: backscatter.dmap
    :members:
    :undoc-members:
    :show-inheritance:
