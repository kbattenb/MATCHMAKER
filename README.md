MATCHMAKER
--------------------
Given a motif sequence(s) and a reference dataset, runs an HMM search iteratively until the search converges.

Please cite:<br/>
(Unpublished at the moment)



DEPENDENCIES
--------------------
(required)<br/>
POSIX (Perl module)<br/>
Math::CDF (Perl module)<br/>

HMMER (tested version is v3.3)<br/>
MAFFT (tested version is v7.4)<br/>
WebLogo (tested version is v3.6)<br/>

(highly recommended)<br/>
miniconda<br/>



INSTALL
--------------------
MATCHMAKER is a perl script and does not require any installation. However, before running, set the proper path for script directory. (At line 40)

Then simply run MATCHMAKER.pl once like the following:

     perl <SOME/PATH/TO/MATCHMAKER.pl>

This will check for each dependency and reports to you know if there is anything missing. If miniconda is installed, it will also give you the command to install the missing
package or tool.



TEST RUN
--------------------
MATCHMAKER comes with a test dataset.

For a test run, try the following command:

     perl MATCHMAKER.pl --element_name TEST\
                        --query <PATH/TO/MATCHMAKER/TEST/QUERY.fas>\
                        --database <PATH/TO/MATCHMAKER/TEST/DATABASE.fas>\
                        --out_dir <PATH/TO/OUTPUT/DIR>

Also to get help, try the following command:

     perl <SOME/PATH/TO/MATCHMAKER.pl> --help
