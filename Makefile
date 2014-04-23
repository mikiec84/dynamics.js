all:
	coffee -c -o lib/ src/
	uglifyjs lib/dynamics.js -m -c > lib/dynamics.min.js
