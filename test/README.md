# Testing
## BATS

To run a more complete test than the command given in the root README.md, use
bats to run the tests in test/bats. Make sure you have conda and bats installed
and then run:

    $ bats test/bats

The tests will build the necessary conda enviroments, and run the workflow a
few different ways.
