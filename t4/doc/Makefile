# ---------------------------------------------------------
# type "make" command in Unix to create report.pdf 
# ---------------------------------------------------------

#list of LibreOffice Draw files
FILE=report
LODFIGS = $(patsubst %.odg,%.pdf,$(wildcard *.odg))

all: $(FILE).pdf
	evince $< &

$(FILE).pdf: $(wildcard *.tex) $(LODFIGS) T2_octave ngspice1 ngspice2
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
	pdflatex  $(FILE).tex
	
T2_octave:
	octave ../mat/t4.m > octave.log
	
ngspice1: ngspice.log1 
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	
	$(eval NGSFIGS3=$(shell grep _FIG3 $< | sed 's/_FIG3//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS3=$(addsuffix .ps, $(NGSFIGS3)))
	$(eval NGSTABS3=$(shell grep _TAB3 $< | sed 's/_TAB3//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach j, $(NGSTABS3), sed -n '/^$j_TAB3/,/^$j_END3/{p;/^$j_END3/q}' $< | grep -v $j_TAB3 | grep -v $j_END3 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $j_tab3.tex;)
	$(foreach j, $(NGSFIGPS3), ps2pdf $j;)
	
	
	

ngspice2: ngspice.log2
	$(eval NGSFIGS2=$(shell grep _FIG2 $< | sed 's/_FIG2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS2=$(addsuffix .ps, $(NGSFIGS2)))
	$(eval NGSTABS2=$(shell grep _TAB2 $< | sed 's/_TAB2//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach l, $(NGSTABS2), sed -n '/^$l_TAB2/,/^$l_END2/{p;/^$l_END2/q}' $< | grep -v $l_TAB2 | grep -v $l_END2 | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $l_tab2.tex;)
	$(foreach l, $(NGSFIGPS2), ps2pdf $l;)
	
	
ngspice.log1: ../sim/t4.net
	ngspice -b $< -o $@
	
ngspice.log2: ../sim/t4out.net
	ngspice -b $< -o $@
	



#convert libreoffice draw figures to pdf
%.pdf: %.odg
	soffice --headless --convert-to pdf $<

clean:
	@rm -rf *.aux *.bbl *.blg *.glg *.glo *.gls *.ilg *.ist *.lof *.txt *.log1 *.log2 *.log3 *.log4 *.log5 *tab2.tex
	@rm -rf *.log *.lot *.nlo *.nls *.out *.toc *~ *.*% ./*.pdf ./*.ps ./*.txt
	@rm -rf *_tab.tex octave-workspace *.eps
	@rm -f ../doc/TA1.tex
	@rm -f ../doc/TA2.tex
	@rm -f ../doc/TA3.tex
	@rm -f ../doc/TA4.tex
	@rm -f ../doc/TA5.tex
	@rm -f ../doc/op_TAB.tex
	@rm -f ../doc/op_TAB2.tex
	@rm -f ../doc/op_tab3.tex

.PHONY: all clean ngspice octave
