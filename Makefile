
all: fig2 fig3 fig5

fig2:
	scripts/fig2.sh ./

fig3:
	scripts/fig3.sh ./

fig5:
	scripts/fig5.sh ./

clean:
	rm -rf output/ temp/

.PHONY: all fig2 fig3 fig5 clean
