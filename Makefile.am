SUBDIRS = libpssa
if LIBSBML
SUBDIRS += libpssa pssa_cli pssa_test examples
endif
SUBDIRS += doc

EXTRA_DIST = 
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory

# include Doxygen rules (requires autoconf-archive >2016-03-20)
@DX_RULES@

clean-local:
	-rm -rf $(DX_DOCDIR)/*

doxygen-prepare:
	$(MKDIR_P) $(DX_DOCDIR)

docs: doxygen-prepare doxygen



# ARTICLE_FILES = 

# install-data-local:
# 	$(MKDIR_P) $(DESTDIR)$(datadir)/$(PACKAGE)/tutorial/media
# 	#for f in $(ARTICLE_FILES); do \
# 	#	$(INSTALL_DATA) $$f \
# 	#	$(DESTDIR)$(datadir)/$(PACKAGE)/; done
# 	$(INSTALL_DATA) ./doc/tutorial/tutorial.html \
# 		$(DESTDIR)$(datadir)/$(PACKAGE)/tutorial
# 	for f in $(TUTORIAL_MEDIA_FILES); do \
# 		$(INSTALL_DATA) $$f \
# 		$(DESTDIR)$(datadir)/$(PACKAGE)/tutorial/media; done
# 
# uninstall-local:
# 	for f in $(TUTORIAL_MEDIA_FILES); do \
# 		rm $(DESTDIR)$(datadir)/$(PACKAGE)/tutorial/media/$$f; done
# 	#for f in $(ARTICLE_FILES); do \
# 	#	rm $(DESTDIR)$(datadir)/$(PACKAGE)/$$f; done
# 	rm $(DESTDIR)$(datadir)/$(PACKAGE)/tutorial/tutorial.html
# 	-rmdir $(DESTDIR)$(datadir)/$(PACKAGE)/tutorial/media
# 	-rmdir $(DESTDIR)$(datadir)/$(PACKAGE)/tutorial


