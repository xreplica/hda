# hda

This is a re-working using better tools and practices of my old implementation of the HDA (for Heuristic Divisive Analysis) clustering algorithm described in http://scitepress.org/DigitalLibrary/Link.aspx?doi=10.5220/0004925902010210


# Requirements
* C++11 compiler (tested with GCC 4.8.2)
* CMake >=2.8.2

# Compiling
Making the project:
`mkdir build`
`cd build`
`cmake ..`
`make`
To run unit tests:
`make test`

# Execution
The executable is nammed "hda". To run the clustering method, run the following command:
$ ./hda [options]
	required flags:
-f <dataset filename>
-d <dataset vector dimension>
	optional flags:
-m <merging mathod> (0=none, 1=nearest neighbor(default), 2=k-means+/chebychev)
-l <merging/linking parameter>
-s <# of iterations to stop splitting>
-D <distance metric> (0=euclidean(default) , 1=mahalanobis(deprecated), 2=gaussian kernel)
-M (no splitting/merge only)
-u (display partition matrix)
-k (display final clusters)
-i (display validation indices)
-v (verbose output)
-h (display help)

# Other clustering algorithms
Along with HDA, multiple other clustering algorithms as described in the article are implemented as well.
