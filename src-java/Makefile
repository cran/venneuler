JAVAC=javac
JAR=jar
JAVADOC=javadoc
JDKVER=1.5
JFLAGS=-source $(JDKVER) -target $(JDKVER)
JAPIURL=http://download.oracle.com/javase/1.5.0/docs/api/

DOCSDIR=javaDoc

SRC=edu/uic/ncdm/venn/data/FilePicker.java edu/uic/ncdm/venn/data/FileReader.java edu/uic/ncdm/venn/data/VennData.java edu/uic/ncdm/venn/display/VennCanvas.java edu/uic/ncdm/venn/display/VennFrame.java edu/uic/ncdm/venn/Eigen.java edu/uic/ncdm/venn/Venn.java edu/uic/ncdm/venn/VennAnalytic.java edu/uic/ncdm/venn/VennDiagram.java

all: venneuler.jar

venneuler.jar: $(SRC)
	@rm -rf obj
	@mkdir obj
	$(JAVAC) -d obj $(JFLAGS) $(SRC)
	@for i in $(SRC); do cp $$i obj/$$i; done
	(cd obj; $(JAR) fc ../$@ edu; cd ..)
	@rm -rf obj

docs: $(SRC)
	rm -rf $(DOCSDIR)
	mkdir $(DOCSDIR)
	$(JAVADOC) -d $(DOCSDIR) -author -version -breakiterator -link $(JAPIURL) $(SRC)

install: venneuler.jar
	-mkdir -p ../inst/java
	cp venneuler.jar ../inst/java/venneuler.jar

clean:
	rm -rf obj edu/uic/ncdm/venn/*.class edu/uic/ncdm/venn/*/*.class venneuler.jar $(DOCSDIR)
