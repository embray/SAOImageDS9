#  Copyright (C) 1999-2021
#  Smithsonian Astrophysical Observatory, Cambridge, MA, USA
#  For conditions of distribution and use, see copyright notice in "copyright"

package provide DS9 1.0

# Native ASDF support (see ASDF_NATIVE_SUPPORT_DESIGN.md, TODO.md).
#
# The container itself - the magic line, the YAML tree, the block layout
# and headers, the codecs, the ndarray enumeration, the datatype mapping
# and the mask resolution - is read in C++, in fitsy/asdf.C, the way every
# other format DS9 reads is. What is left here is the part that
# is genuinely text manipulation: lifting a GWCS subtree out of the tree
# and reshaping it into a standalone document AST will accept, and the
# UI/scripting plumbing around it.
#
# Two doors into that C++ reader:
#
#   `$frame load asdf {name} {file} {path} layer' - the pixel path. The
#   frame builds a FitsAsdf (fitsy/asdf.C) which seeks to the array's
#   block, decompresses it, widens a datatype fitsy cannot hold, resolves
#   an asdf mask onto FITS BLANK/NaN, and hands the frame a synthetic
#   header. No temp file, no Tcl byte array, and the file is never read
#   whole - which is the point, since a Roman *_cal.asdf is ~200MB and
#   used to be slurped into a Tcl string in its entirety.
#
#   `fitsy asdf tree|arrays|block <file>' - the metadata path, through
#   tclfitsy. Same C++ reader, so there is one implementation of the
#   format rather than one for pixels and another for everything else.

# Every block backed core/ndarray in `fn', as a list of
#   {path source datatype byteorder shapelist unsupported maskscalar bitpix}
# `path' is root relative (roman/data,
# roman/meta/wcs/steps/0/transform/forward/1/coefficients); `unsupported'
# names the reasons the array cannot be read as pixels, empty when there
# are none; `bitpix' is what the loader would use after any lossless
# widening, empty for a datatype it cannot represent.
#
# Raises if `fn' is not an ASDF file or cannot be read.
proc AsdfArrays {fn} {
    return [fitsy asdf arrays $fn]
}

# The file's YAML tree text: everything before the first binary block.
proc AsdfTree {fn} {
    return [fitsy asdf tree $fn]
}

# Looks up one enumerated node by its root-relative path in an
# AsdfArrays result. Returns the entry minus its path -
# {source datatype byteorder shapelist unsupported maskscalar bitpix} -
# or {} if there is no such path.
proc AsdfFindNdarrayPath {entries path} {
    foreach entry $entries {
	if {[lindex $entry 0] eq $path} {
	    return [lrange $entry 1 end]
	}
    }
    return {}
}

# Turns a user-supplied key into a full path. A key containing "/" is
# already one. A bare name means Roman's fixed top-level layout
# (roman/data, roman/err, ...); failing that, a unique match on the last
# path component is accepted, so a non-Roman file's `sci` resolves
# without the user having to spell out its full nesting. An ambiguous
# bare name returns {} rather than picking one.
proc AsdfResolvePath {entries key} {
    if {[string first / $key] >= 0} {
	return $key
    }

    foreach entry $entries {
	if {[lindex $entry 0] eq "roman/$key"} {
	    return "roman/$key"
	}
    }

    set hits {}
    foreach entry $entries {
	set pp [lindex $entry 0]
	if {[file tail $pp] eq $key} {
	    lappend hits $pp
	}
    }
    if {[llength $hits] == 1} {
	return [lindex $hits 0]
    }

    return {}
}

# Native ASDF/GWCS WCS attachment (TODO.md Phase 3). Ported from the Phase 1
# Python spike (utils/asdf_gwcs_probe/extract_subtree.py,
# resolve_ndarray.py) - see PHASE3_WCS_HANDOFF.md for the full design and
# why this goes in via the frame's existing `wcs replace <which> <filename>`
# command (a temp file, NOT `wcs replace text ...` - the latter's STRING
# token is unescaped and truncates on YAML's own literal "/{/} characters).

# Finds `$indent$key: ...` (default 2-space indent - a direct child of the
# top-level "roman:" mapping) and returns its dedented, standalone subtree
# text (still tagged with its original key name, not yet renamed to
# "wcs:"). Returns {} if not found. Mirrors extract_subtree.py's
# extract_key_subtree. `indent` is a parameter (not hardcoded to 2 spaces)
# because real Roman files put the WCS at different nesting depths
# depending on product type - see AsdfExtractWcsText below.
proc AsdfExtractKeySubtree {tree key {indent "  "}} {
    set marker "\n$indent$key:"
    set start [string first $marker $tree]
    if {$start < 0} {
	return {}
    }
    set bodyStart [expr {$start + 1}]
    set bodyEnd [string length $tree]
    set nextKeyRe "\n$indent\[A-Za-z_\]\[A-Za-z0-9_\]*:"
    if {[regexp -indices -start $bodyStart -- $nextKeyRe $tree m]} {
	set bodyEnd [lindex $m 0]
    }
    set subtree [string trimright \
		     [string range $tree $bodyStart [expr {$bodyEnd - 1}]]]
    if {[string range $subtree end-2 end] eq "..."} {
	set subtree [string trimright [string range $subtree 0 end-3]]
    }
    set n [string length $indent]
    set lines {}
    foreach line [split $subtree "\n"] {
	if {[string range $line 0 [expr {$n - 1}]] eq $indent} {
	    lappend lines [string range $line $n end]
	} else {
	    lappend lines $line
	}
    }
    return "[join $lines "\n"]\n"
}

# A proper GWCS object's tag is `!<tag:...:gwcs/wcs-...>`; a bare transform
# (e.g. a pre-launch/pysiaf distortion file's coordinate_distortion_transform,
# tagged !transform/...) has no such wrapper and needs one synthesized -
# astRead() rejects a bare transform as a document root outright (TODO.md
# Phase 1).
proc AsdfIsBareTransform {subtree} {
    set firstLine [lindex [split $subtree "\n"] 0]
    return [expr {![string match "*gwcs/wcs*" $firstLine]}]
}

# `key: !<gwcs/wcs-...> ...` -> `wcs: !<gwcs/wcs-...> ...` (rename the
# top-level key only; astRead() dispatches on this shape directly).
proc AsdfAsWcsShell {key subtree} {
    set idx [string first "$key:" $subtree]
    if {$idx < 0} {
	return $subtree
    }
    return [string replace $subtree $idx [expr {$idx + [string length $key]}] \
		"wcs:"]
}

# Wraps a bare `key: !transform/...` value in a minimal synthetic 2-step
# GWCS shell (detector frame carrying the real transform, feeding a bare
# terminal frame) that astRead() does accept. Mirrors
# extract_subtree.py's SHELL_TEMPLATE/as_wrapped_transform exactly.
proc AsdfWrapBareTransform {subtree} {
    set colon [string first ":" $subtree]
    set valueOnly [string trimleft [string range $subtree \
					 [expr {$colon + 1}] end]]
    set lines [split [string trimright $valueOnly "\n"] "\n"]
    set rest {}
    foreach line [lrange $lines 1 end] {
	if {[string trim $line] eq {}} {
	    lappend rest $line
	} else {
	    lappend rest "    $line"
	}
    }
    set reindented "[lindex $lines 0]\n[join $rest "\n"]"

    set template {wcs: !<tag:stsci.edu:gwcs/wcs-1.4.0>
  name: distortion_shell
  steps:
  - !<tag:stsci.edu:gwcs/step-1.3.0>
    frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>
      axes_names: [x, y]
      axes_order: [0, 1]
      axis_physical_types: ["custom:x", "custom:y"]
      name: detector
      unit: [!unit/unit-1.0.0 pixel, !unit/unit-1.0.0 pixel]
    transform: @TRANSFORM@
  - !<tag:stsci.edu:gwcs/step-1.3.0>
    frame: !<tag:stsci.edu:gwcs/frame2d-1.2.0>
      axes_names: [v2, v3]
      axes_order: [0, 1]
      axis_physical_types: ["custom:v2", "custom:v3"]
      name: output
      unit: [!unit/unit-1.0.0 arcsec, !unit/unit-1.0.0 arcsec]
}
    return [string map [list @TRANSFORM@ $reindented] $template]
}

# ASDF/numpy datatype+byteorder -> a Tcl `binary scan` format spec (see
# tcl9.0/doc/binary.n: lowercase little-endian/uppercase big-endian, `u`
# suffix unsigned, `r`/`R` float32, `q`/`Q` float64). Returns {} for
# unsupported types. Used only on the small inline coefficient arrays a
# GWCS subtree carries, never on pixels.
proc AsdfBinaryFmt {datatype byteorder count} {
    set big [expr {$byteorder eq "big"}]
    switch -- $datatype {
	float64 {set code [expr {$big ? "Q" : "q"}]}
	float32 {set code [expr {$big ? "R" : "r"}]}
	int64   {set code [expr {$big ? "W" : "w"}]}
	uint64  {set code [expr {$big ? "Wu" : "wu"}]}
	int32   {set code [expr {$big ? "I" : "i"}]}
	uint32  {set code [expr {$big ? "Iu" : "iu"}]}
	int16   {set code [expr {$big ? "S" : "s"}]}
	uint16  {set code [expr {$big ? "Su" : "su"}]}
	int8    {set code "c"}
	uint8   {set code "cu"}
	default {return {}}
    }
    return "$code$count"
}

# Formats a flat list of decoded numeric values as a YAML-nested list
# matching `shape`, indented under `indent`. Mirrors resolve_ndarray.py's
# _nest.
proc AsdfNest {vals shape indent} {
    if {[llength $shape] <= 1} {
	return "\[[join $vals {, }]\]"
    }
    set cols 1
    foreach s [lrange $shape 1 end] {
	set cols [expr {$cols * $s}]
    }
    set lines {}
    for {set r 0} {$r < [lindex $shape 0]} {incr r} {
	set rowVals [lrange $vals [expr {$r*$cols}] [expr {$r*$cols+$cols-1}]]
	if {[llength $shape] == 2} {
	    lappend lines "$indent- \[[join $rowVals {, }]\]"
	} else {
	    lappend lines "$indent- [AsdfNest $rowVals [lrange $shape 1 end] \
					  "$indent  "]"
	}
    }
    return "\n[join $lines "\n"]"
}

# Decodes `decoded` (raw bytes for one already-decompressed block) per
# `datatype`/`byteorder`/`shapelist` and nests it as a YAML `data:` literal.
proc AsdfNestValues {decoded datatype byteorder shapelist indent} {
    set n 1
    foreach s $shapelist {
	set n [expr {$n * $s}]
    }
    set fmt [AsdfBinaryFmt $datatype $byteorder $n]
    if {$fmt eq {}} {
	error "unsupported ndarray datatype: $datatype"
    }
    binary scan $decoded $fmt vals
    return [AsdfNest $vals $shapelist $indent]
}

# Rewrites every block-sourced `!core/ndarray-*` node in `text` (a GWCS
# subtree, already extracted/renamed/wrapped) as an inline `data:` literal,
# reading each one's real block bytes out of `fn` as it goes - AST's
# YamlChan only accepts inline ndarrays, never block-sourced ones (design
# doc SS7b). Real Roman GWCS documents always use block-sourced ndarrays,
# even for tiny coefficient matrices - this is required, not an edge case
# (TODO.md Phase 1). Mirrors resolve_ndarray.py's resolve().
proc AsdfResolveNdarrays {text fn} {
    # Two Tcl-regex-specific gotchas fixed here, both confirmed the hard
    # way against a real multi-hundred-line GWCS subtree, not assumed:
    # (1) Tcl's "." matches newlines too, unlike most other regex flavors
    #     (no -line option is in play here) - a naive ".*\n" per mid-line
    #     doesn't stop at end-of-line at all, letting a single "iteration"
    #     swallow arbitrarily many real lines. Fixed by using "[^\n]*"
    #     instead of ".*" so each iteration is truly one physical line.
    # (2) Tcl's ARE engine picks an overall leftmost-LONGEST match (POSIX
    #     semantics), not Perl-style leftmost-shortest, so even a properly
    #     line-bounded but unrestricted "*" mid-loop would still prefer
    #     matching through many more (non-"shape:") lines than intended
    #     before hitting a later, wrong "shape:" occurrence. Fixed with an
    #     explicit "(?!shape:)" negative lookahead so a "shape:" line can
    #     never be consumed by the generic loop in the first place, which
    #     removes the ambiguity outright instead of relying on greediness.
    set re {!core/ndarray-([0-9.]+)\n([ \t]+)source:[ \t]*([0-9]+)[ \t]*\n((?:[ \t]+(?!shape:)[A-Za-z_]+:[^\n]*\n)*)([ \t]+)shape:[ \t]*\[([0-9, ]+)\][ \t]*\n}
    set result {}
    set pos 0
    while {[regexp -indices -start $pos -- $re $text m verIdx indentIdx \
		srcIdx midIdx indent2Idx shapeIdx]} {
	lassign $m mStart mEnd
	append result [string range $text $pos [expr {$mStart - 1}]]

	set ver [string range $text {*}$verIdx]
	set indent [string range $text {*}$indentIdx]
	set source [string range $text {*}$srcIdx]
	set mid [string range $text {*}$midIdx]
	set indent2 [string range $text {*}$indent2Idx]
	set shapeText [string range $text {*}$shapeIdx]

	set datatype float64
	regexp {datatype:[ \t]*([A-Za-z0-9_]+)} $mid -> datatype
	set byteorder little
	regexp {byteorder:[ \t]*([A-Za-z0-9_]+)} $mid -> byteorder

	set shapelist {}
	foreach s [split $shapeText ,] {
	    lappend shapelist [string trim $s]
	}

	set decoded [fitsy asdf block $fn $source]
	set dataRepr [AsdfNestValues $decoded $datatype $byteorder \
			  $shapelist $indent]

	append result "!core/ndarray-$ver\n${indent}data: $dataRepr\n$mid${indent2}shape: \[$shapeText\]\n"
	set pos [expr {$mEnd + 1}]
    }
    append result [string range $text $pos end]
    return $result
}

# Top-level Phase 3 entry point: extracts a standalone, AST-ready `wcs:`
# YAML document from the file's tree text, trying each candidate key in
# turn (the key name/nesting isn't fixed across Roman product types -
# TODO.md Phase 1). Returns {} if none of the candidate keys are found.
# `fn` is the ASDF file, needed to resolve the subtree's own embedded
# ndarrays back into inline literals.
#
# Two real, confirmed shapes, tried in order:
#  - `keys`, at 2-space indent directly under "roman:" - the small
#    WCS-only distribution products (e.g. Roman Data Workshop's
#    "*_wcs.asdf" files) use this for wcs_l1/wcs_l2.
#  - `metaKeys`, at 4-space indent nested under "roman: {meta: {...}}" -
#    the design doc's original roman.meta.wcs assumption, confirmed
#    against a real full science-product "*_cal.asdf" file.
proc AsdfExtractWcsText {tree fn {keys {wcs wcs_l2 wcs_l1}} {metaKeys wcs}} {
    # Indents to try, in order. Two-space first because that is where every
    # real Roman product keeps it (`roman:` then `  meta:` ... or `  wcs:`),
    # so the common case still matches on the first attempt.
    #
    # The empty indent - a `wcs:` key at the very top of a flat tree - was
    # missing until fixtures for the sky projections turned up:
    # Phase 4 generalized the *array* lookup away from Roman's fixed
    # indents (the ndarray enumerator walks the whole tree) but left the *WCS*
    # lookup pinned to them, so a perfectly legal flat GWCS file loaded its
    # pixels and silently got no WCS at all.
    set subtree {}
    foreach indent {"  " ""} {
	foreach key $keys {
	    set subtree [AsdfExtractKeySubtree $tree $key $indent]
	    if {$subtree ne {}} {
		break
	    }
	}
	if {$subtree ne {}} {
	    break
	}
    }
    if {$subtree eq {}} {
	foreach key $metaKeys {
	    set subtree [AsdfExtractKeySubtree $tree $key "    "]
	    if {$subtree ne {}} {
		break
	    }
	}
    }
    if {$subtree eq {}} {
	return {}
    }

    if {[AsdfIsBareTransform $subtree]} {
	set doc [AsdfWrapBareTransform $subtree]
    } else {
	set doc [AsdfAsWcsShell $key $subtree]
    }

    # The standalone document needs its own "%TAG !" directive - real GWCS
    # subtrees use the short-form tag shorthand (e.g. !transform/compose-
    # 1.4.0) throughout, not just at the top level, and that shorthand only
    # expands to the full tag:stsci.edu:asdf:... form in the presence of
    # this directive. Mirrors extract_subtree.py's HEADER constant exactly.
    set header "#ASDF 1.0.0\n%YAML 1.1\n%TAG ! tag:stsci.edu:asdf/\n--- !core/asdf-1.1.0\n"

    return [AsdfResolveNdarrays "$header$doc" $fn]
}

# --------------------------------------------------------------------
# Header/metadata display (TODO.md Phase 3, design doc SS8 point 4).
#
# An ASDF file's YAML tree *is* its header, and it is already exactly the
# right thing to show verbatim: ASDF never writes an array inline, only a
# `source: N` reference alongside its datatype/shape, so even the 197MB
# Roman cal file's tree is ~60KB of readable YAML. So there is no card
# list to build and no "report arrays by shape/dtype" pass to write - the
# format did that for us.
#
# Cached per frame rather than re-read on demand: the C++ reader can hand
# back a tree any time (`fitsy asdf tree'), but the *frame* is what the
# viewer is asked about, and a frame remembers which file it came from
# only through this cache. Cleared whenever that frame is reloaded or
# unloaded (see AsdfClearTree's callers in load.tcl and frame.tcl), so a
# FITS file loaded over an ASDF one can never show the previous file's
# tree. The frame's own header is synthetic - the dimensions and BITPIX
# FitsAsdf gave it, plus BLANK - and has nothing worth displaying.

proc AsdfSetTree {frame fn tree} {
    global asdf

    if {$frame == {}} {
	return
    }
    set asdf(tree,$frame) $tree
    set asdf(file,$frame) $fn
}

proc AsdfHasTree {frame} {
    global asdf

    return [expr {$frame != {} && [info exists asdf(tree,$frame)]}]
}

proc AsdfGetTree {frame} {
    global asdf

    if {![AsdfHasTree $frame]} {
	return {}
    }
    return $asdf(tree,$frame)
}

# Closes the viewer window without dropping the cache - what the `header
# close` command wants. DestroyHeader/DestroyHeaderOne in header.tcl only
# know about the hd-* FITS header windows, so this is a separate teardown.
proc AsdfDestroyHeader {frame} {
    if {$frame == {}} {
	return
    }

    set varname "asdf-$frame"
    global $varname
    if {[info exists $varname]} {
	SimpleTextDestroy $varname
    }

    AsdfHeaderTreeDestroy "asdfhdr-$frame"
}

# Drops the cache and closes the viewer with it, so it can't keep showing
# a file the frame no longer holds.
proc AsdfClearTree {frame} {
    global asdf

    if {$frame == {}} {
	return
    }
    unset -nocomplain asdf(tree,$frame)
    unset -nocomplain asdf(file,$frame)

    AsdfDestroyHeader $frame
}

# Shows the cached tree in the same SimpleTextDialog widget the FITS
# header viewer uses (header.tcl's DisplayHeader), so it inherits that
# window's Save/Print/Find machinery for free. No keyword tagging: that
# loop tags a fixed 8-character column, which is a FITS card convention
# with no YAML equivalent.
proc DisplayAsdfHeader {frame} {
    global asdf

    if {![AsdfHasTree $frame]} {
	return
    }

    # the tree navigator; the plain text dump if the YAML cannot be walked
    if {[catch {fitsy asdf yaml $asdf(file,$frame)} nodes]} {
	AsdfHeaderText $frame
	return
    }
    AsdfHeaderTreeDialog $frame $nodes
}

# The YAML tree as plain text, which is also what the header viewer
# showed before it became a tree.
proc AsdfHeaderText {frame} {
    global asdf

    if {![AsdfHasTree $frame]} {
	return
    }

    set varname "asdf-$frame"
    global $varname

    SimpleTextDialog $varname [file tail $asdf(file,$frame)] 80 40 \
	insert top $asdf(tree,$frame)
}

# The header viewer for an ASDF frame: the YAML directives as text, then
# the document as a tree. Any node can be copied, which copies its source
# text verbatim (see AsdfYamlNodeText). `nodes' is `fitsy asdf yaml'.
proc AsdfHeaderTreeDialog {frame nodes} {
    global asdf
    global asdftree
    global ds9

    set varname "asdfhdr-$frame"
    upvar #0 $varname var
    global $varname

    set var(top) ".${varname}"
    set var(mb) ".${varname}mb"
    if {[winfo exists $var(top)]} {
	raise $var(top)
	return
    }

    set var(frame) $frame
    set tree $asdf(tree,$frame)
    set var(tree) $tree

    set w $var(top)
    set mb $var(mb)
    Toplevel $w $mb 7 [file tail $asdf(file,$frame)] \
	[list AsdfHeaderTreeDestroy $varname]

    # Save and Print are SimpleText's, writing the tree text rather than a
    # text widget's contents, and saving as YAML
    set var(contents) AsdfHeaderTreeContents
    set var(fbox) yamlfbox
    set var(search) {}

    $mb add cascade -label [msgcat::mc {File}] -menu $mb.file
    ThemeMenu $mb.file
    $mb.file add command -label "[msgcat::mc {Save}]..." \
	-command [list AsdfHeaderTreeSave $varname] \
	-accelerator "${ds9(ctrl)}S"
    $mb.file add separator
    $mb.file add command -label "[msgcat::mc {Print}]..." \
	-command [list SimpleTextPrint $varname] \
	-accelerator "${ds9(ctrl)}P"
    $mb.file add separator
    $mb.file add command -label "[msgcat::mc {View as Text}]..." \
	-command [list AsdfHeaderText $frame]
    $mb.file add separator
    $mb.file add command -label [msgcat::mc {Close}] \
	-command [list AsdfHeaderTreeDestroy $varname] \
	-accelerator "${ds9(ctrl)}W"

    $mb add cascade -label [msgcat::mc {Edit}] -menu $mb.edit
    ThemeMenu $mb.edit
    $mb.edit add command -label [msgcat::mc {Copy}] \
	-command [list AsdfHeaderTreeCopy $varname] \
	-accelerator "${ds9(ctrl)}C"
    $mb.edit add separator
    $mb.edit add command -label "[msgcat::mc {Find}]..." \
	-command [list AsdfHeaderTreeFind $varname] \
	-accelerator "${ds9(ctrl)}F"
    $mb.edit add command -label [msgcat::mc {Find Next}] \
	-command [list AsdfHeaderTreeFindNext $varname] \
	-accelerator "${ds9(ctrl)}G"
    $mb.edit add separator
    $mb.edit add command -label [msgcat::mc {Expand All}] \
	-command [list AsdfHeaderTreeOpen $varname 1]
    $mb.edit add command -label [msgcat::mc {Collapse All}] \
	-command [list AsdfHeaderTreeOpen $varname 0]

    # the directives: everything before the document's `---'
    set directives {}
    foreach line [split $tree \n] {
	if {[string match ---* $line]} {
	    break
	}
	lappend directives $line
    }
    set var(text) $w.text
    roText::roText $var(text)
    $var(text) configure -height [llength $directives] -width 80 \
	-wrap none -fg [ThemeTreeForeground] -bg [ThemeTreeBackground] \
	-state normal
    $var(text) insert end [join $directives \n]

    # the document
    set f [ttk::frame $w.tf]
    set tv [ttk::treeview $f.tree \
		-yscroll [list $f.yscroll set] \
		-xscroll [list $f.xscroll set] \
		-selectmode extended \
		-columns {tag value} \
		-height 24 \
	       ]
    set var(tv) $tv
    ttk::scrollbar $f.yscroll -command [list $tv yview]
    ttk::scrollbar $f.xscroll -command [list $tv xview] -orient horizontal
    $tv heading #0 -text [msgcat::mc {Node}] -anchor w
    $tv heading tag -text [msgcat::mc {Tag}] -anchor w
    $tv heading value -text [msgcat::mc {Value}] -anchor w
    $tv column #0 -width 260 -stretch 0
    $tv column tag -width 260 -stretch 0
    $tv column value -width 360 -stretch 1

    AsdfFillYamlTree $varname $nodes

    grid $tv $f.yscroll -sticky news
    grid $f.xscroll -sticky ew
    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 0 -weight 1

    pack $var(text) -side top -fill x
    pack [ttk::separator $w.sep -orient horizontal] -side top -fill x
    pack $f -side top -fill both -expand true

    bind $tv <<Copy>> [list AsdfHeaderTreeCopy $varname]
    bind $w <<Close>> [list AsdfHeaderTreeDestroy $varname]
    bind $w <<Save>> [list AsdfHeaderTreeSave $varname]
    bind $w <<Print>> [list SimpleTextPrint $varname]
    bind $w <<Find>> [list AsdfHeaderTreeFind $varname]
    bind $w <<FindNext>> [list AsdfHeaderTreeFindNext $varname]
}

proc AsdfHeaderTreeContents {varname} {
    upvar #0 $varname var
    global $varname

    return $var(tree)
}

proc AsdfHeaderTreeSave {varname} {
    upvar #0 $varname var
    global $varname
    global asdf
    global yamlfbox

    # suggest the ASDF file's own name, as YAML
    set yamlfbox(file) "[file rootname [file tail $asdf(file,$var(frame))]].yaml"
    SimpleTextSaveFile $varname
}

# Find searches the nodes in document order (depth first), matching a
# node's key (in a mapping), its tag (only if explicit) and a scalar's
# value, as a case-insensitive regexp like the text view's Find.
proc AsdfHeaderTreeFind {varname} {
    upvar #0 $varname var
    global $varname

    set result $var(search)
    if {[EntryDialog [msgcat::mc {Search}] \
	     [msgcat::mc {Enter Search Expression}] 40 result]} {
	# we might have been closed while waiting
	if {![info exists $varname]} {
	    return
	}
	set var(search) $result
	AsdfHeaderTreeSearch $varname 0
    }
}

proc AsdfHeaderTreeFindNext {varname} {
    upvar #0 $varname var
    global $varname

    if {$var(search) eq {}} {
	return
    }
    set sel [lindex [$var(tv) selection] 0]
    AsdfHeaderTreeSearch $varname [expr {$sel eq {} ? 0 : $sel + 1}]
}

# Selects the first node at or after `from' that matches, wrapping around.
proc AsdfHeaderTreeSearch {varname from} {
    upvar #0 $varname var
    global $varname

    set nn $var(nnodes)
    for {set ii 0} {$ii < $nn} {incr ii} {
	set id [expr {($from + $ii) % $nn}]
	foreach str $var(match,$id) {
	    if {[catch {regexp -nocase -- $var(search) $str} hit]} {
		Error "[msgcat::mc {Invalid search expression}]: $var(search)"
		return
	    }
	    if {$hit} {
		$var(tv) see $id
		$var(tv) selection set [list $id]
		$var(tv) focus $id
		return
	    }
	}
    }
    Error "$var(search) [msgcat::mc {Not Found}]"
}

proc AsdfHeaderTreeDestroy {varname} {
    upvar #0 $varname var
    global $varname

    if {[info exists var(top)] && [winfo exists $var(top)]} {
	destroy $var(top)
	destroy $var(mb)
    }
    unset -nocomplain $varname
}

# Inserts every node of `fitsy asdf yaml' into the dialog's treeview, item
# id being the node's index. Collections show their size, scalars their
# value on one line; the full text is what Copy gives.
proc AsdfFillYamlTree {varname nodes} {
    upvar #0 $varname var
    global $varname
    global asdftree

    set tv $var(tv)

    set nchild {}
    foreach node $nodes {
	dict incr nchild [lindex $node 0]
    }

    set ii 0
    foreach node $nodes {
	lassign $node parent kind label tag anchor value keystart start end

	if {$parent < 0} {
	    set pid {}
	    set depth 0
	    set label [msgcat::mc {(document)}]
	} else {
	    set pid $parent
	    set depth [expr {$var(depth,$parent) + 1}]
	}
	set var(depth,$ii) $depth
	set var(span,$ii) [list $keystart $start $end]

	# what Find matches: the key, an explicit tag, a scalar's value
	set var(match,$ii) {}
	if {$parent >= 0 && [lindex $nodes $parent 1] eq {m}} {
	    lappend var(match,$ii) $label
	}
	if {$tag ne {}} {
	    lappend var(match,$ii) $tag
	}
	if {$kind eq {v}} {
	    lappend var(match,$ii) $value
	}

	if {[lindex $nodes $parent 1] eq {s}} {
	    set label "\[$label\]"
	}

	set nn [expr {[dict exists $nchild $ii] ? [dict get $nchild $ii] : 0}]
	switch -- $kind {
	    m {set vv "{$nn}"}
	    s {set vv "\[$nn\]"}
	    a {set vv "*$value"}
	    default {
		set vv [string map [list \n { } \r {} \t { }] $value]
		if {[string length $vv] > 120} {
		    set vv "[string range $vv 0 116]..."
		}
	    }
	}
	if {$anchor ne {}} {
	    append tag " &$anchor"
	}

	$tv insert $pid end -id $ii -text $label -values [list $tag $vv] \
	    -open [expr {$depth < $asdftree(headerdepth)}]
	incr ii
    }
    set var(nnodes) $ii
}

proc AsdfHeaderTreeOpen {varname open} {
    upvar #0 $varname var
    global $varname

    set todo [$var(tv) children {}]
    while {$todo ne {}} {
	set todo [lassign $todo id]
	$var(tv) item $id -open $open
	lappend todo {*}[$var(tv) children $id]
    }
}

# Copies the selected nodes' source text; with nothing selected in the tree
# but a selection in the directives, copies that instead.
proc AsdfHeaderTreeCopy {varname} {
    upvar #0 $varname var
    global $varname

    set sel [$var(tv) selection]
    if {$sel eq {}} {
	if {[$var(text) tag ranges sel] ne {}} {
	    clipboard clear -displayof $var(top)
	    clipboard append -displayof $var(top) \
		[$var(text) get sel.first sel.last]
	}
	return
    }

    set out {}
    foreach id $sel {
	lappend out [AsdfYamlNodeText $var(tree) $id {*}$var(span,$id)]
    }
    clipboard clear -displayof $var(top)
    clipboard append -displayof $var(top) [join $out {}]
}

# The source text of one node, from its key (if it has one) to its end,
# de-indented so that it stands on its own as YAML. Node 0, the document,
# is the whole tree.
proc AsdfYamlNodeText {tree id keystart start end} {
    if {$id == 0} {
	return $tree
    }

    set from [expr {$keystart >= 0 ? $keystart : $start}]
    set ls [expr {[string last \n $tree [expr {$from - 1}]] + 1}]
    set col [expr {$from - $ls}]

    # indentation only before the node on its line: take the whole line
    # and strip that indentation from every line. Otherwise (the node
    # follows a `- ' or a flow delimiter) keep the first line as it is.
    set whole [expr {[string trim [string range $tree $ls [expr {$from - 1}]]] eq {}}]
    if {$whole} {
	set from $ls
    }

    set lines [split [string trimright \
			  [string range $tree $from [expr {$end - 1}]]] \n]
    set first {}
    if {!$whole} {
	set first [lindex $lines 0]
	set lines [lrange $lines 1 end]
	set ind $col
	foreach line $lines {
	    if {[string trim $line] eq {}} {
		continue
	    }
	    set nn [expr {[string length $line] - \
			      [string length [string trimleft $line { }]]}]
	    if {$nn < $ind} {
		set ind $nn
	    }
	}
	# A key that starts mid-line sits in a flow mapping, and what wraps
	# onto further lines is its (flow or plain scalar) value - which,
	# now under a key at column 0, must stay indented.
	if {$keystart >= 0} {
	    set ind [expr {max(0, $ind - 2)}]
	}
	set col $ind
    }

    set out {}
    if {!$whole} {
	lappend out $first
    }
    foreach line $lines {
	set nn [expr {[string length $line] - \
			  [string length [string trimleft $line { }]]}]
	lappend out [string range $line [expr {min($nn, $col)}] end]
    }
    # ends with a newline, like any YAML block - and a literal or folded
    # scalar's final newline is part of its value
    return "[join $out \n]\n"
}

# Hands an extracted, self-contained GWCS document to the current frame
# through the existing `wcs replace <which> <filename>` command. The YAML
# travels via a temp file, deliberately not through `wcs replace text ...`:
# that overload re-lexes its argument with tksao/frame/lex.L, whose STRING
# token rules ("[^"]*" and {[^}]*}) have no escape mechanism and truncate on
# the literal " { } characters every real GWCS document contains. Only the
# temp file's name - short and delimiter-free - reaches the lexer this way.
# See PHASE3_WCS_HANDOFF.md.
#
# `which` is 1-indexed: Base::findAllFits(int which) loops on `while (ptr &&
# which)`, so 0 means "stop immediately" and always returns NULL, which
# surfaces as an empty-message TCL_ERROR. 1 is also ds9/library/wcs.tcl's
# own dwcs(ext) convention.
#
# Errors propagate to the caller rather than being swallowed here - the
# caller decides whether a WCS failure is fatal.
# Which pixel grid does the file's GWCS actually describe?
#
# A GWCS sits in the tree beside the science array it maps, and says nothing
# about its own domain that we can read: R6 in WCS_TEST_PLAN.md - AST drops
# `bounding_box` entirely, so the transform we get back carries no hint of
# the pixel extent it is valid over. So the grid has to come from the tree,
# and the science array is the obvious source. Reuse AsdfResolvePath's
# existing convention for what "data" means (roman/data first, else any
# uniquely matching */data) rather than inventing a second rule.
#
# Returns the science array's shape, or {} if there is no resolvable one.
# --------------------------------------------------------------------
# gwcs/fitswcs_imaging -> FITS WCS cards (WCS_TEST_PLAN J-4).
#
# Build22 L3 coadd products express their WCS as a single
# `!<tag:stsci.edu:gwcs/fitswcs_imaging-1.0.0>` node holding crpix, crval,
# cdelt, pc and a projection, instead of the explicit
# compose/shift/polynomial/gnomonic chain the L2 *_cal.asdf files use.
# yamlchan.c has no handler for that tag, so AST returns no FrameSet and the
# coadd loads with no WCS at all.
#
# But that node *is* a FITS WCS, spelled differently - so rather than
# teaching AST a new tag, translate it into ordinary FITS cards and let DS9's
# own FITS WCS path handle it. There is no approximation anywhere in this:
# unlike the H-7 save-as-FITS problem, this transform really is a TAN.

# The indent of a line, in spaces.
proc AsdfLineIndent {line} {
    regexp {^( *)} $line -> pad
    return [string length $pad]
}

# Every number in a string, in order. Used only on `data:` payload lines,
# which hold nothing else.
proc AsdfNumberList {str} {
    return [regexp -all -inline -- \
		{[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?} $str]
}

# The body of the `gwcs/fitswcs_imaging` node in `yamltext`, or {}.
proc AsdfFitsWcsImagingBody {yamltext} {
    set lines [split $yamltext "\n"]
    set n [llength $lines]
    set start -1
    for {set i 0} {$i < $n} {incr i} {
	if {[string match {*gwcs/fitswcs_imaging*} [lindex $lines $i]]} {
	    set start $i
	    break
	}
    }
    if {$start < 0} {
	return {}
    }
    set base [AsdfLineIndent [lindex $lines $start]]
    set out {}
    for {set j [expr {$start + 1}]} {$j < $n} {incr j} {
	set line [lindex $lines $j]
	if {[string trim $line] eq {}} {
	    continue
	}
	if {[AsdfLineIndent $line] <= $base} {
	    break
	}
	lappend out $line
    }
    return [join $out "\n"]
}

# The inline ndarray payload belonging to `key` within a node body, as a
# flat row-major list of numbers.
#
# asdf writes a small inline ndarray two different ways depending on rank -
# `data: [a, b]` for a vector, but a bare `data:` followed by `- [a, b]`
# rows for a matrix - so both shapes are handled. Returns {} if the key is
# absent, or if its ndarray is block-backed (`source:`) rather than inline,
# which is why the scan stops at the next key at or above `key`'s own indent
# instead of running on and picking up a later key's data.
proc AsdfInlineNdarray {body key} {
    set lines [split $body "\n"]
    set n [llength $lines]
    set ki -1
    for {set i 0} {$i < $n} {incr i} {
	if {[regexp "^ *$key:(\[ \]|\$)" [lindex $lines $i]]} {
	    set ki $i
	    break
	}
    }
    if {$ki < 0} {
	return {}
    }
    set kind [AsdfLineIndent [lindex $lines $ki]]
    for {set j [expr {$ki + 1}]} {$j < $n} {incr j} {
	set line [lindex $lines $j]
	if {[string trim $line] eq {}} {
	    continue
	}
	if {[AsdfLineIndent $line] <= $kind} {
	    break
	}
	if {![regexp {^ *data:(.*)$} $line -> rest]} {
	    continue
	}
	if {[string trim $rest] ne {}} {
	    return [AsdfNumberList $rest]
	}
	set out {}
	for {set k [expr {$j + 1}]} {$k < $n} {incr k} {
	    set row [lindex $lines $k]
	    if {![regexp {^ *- *\[} $row]} {
		break
	    }
	    foreach v [AsdfNumberList $row] {
		lappend out $v
	    }
	}
	return $out
    }
    return {}
}

# RADESYS/CTYPE pair for the document's celestial frame. Only the
# equatorial frames are translated: galactic and ecliptic would need
# GLON/GLAT or ELON/ELAT CTYPEs, and emitting an equatorial CTYPE for one of
# those would be silently wrong, so they fall through to {} and the caller
# declines to synthesize.
proc AsdfCelestialRadesys {yamltext} {
    if {![regexp {coordinates/frames/([a-z0-9]+)-} $yamltext -> frame]} {
	return {}
    }
    switch -- $frame {
	icrs {return ICRS}
	fk5 {return FK5}
	fk4 {return FK4}
	default {return {}}
    }
}

# FITS WCS cards equivalent to the document's fitswcs_imaging node, or {}
# if it has none or carries anything this does not translate exactly.
proc AsdfFitsWcsImagingCards {yamltext} {
    set body [AsdfFitsWcsImagingBody $yamltext]
    if {$body eq {}} {
	return {}
    }

    # Only gnomonic (TAN). Every other projection would need its own CTYPE
    # code, and guessing one untested is worse than declining.
    if {![regexp {projection: *!transform/([a-z_0-9]+)-} $body -> proj]} {
	return {}
    }
    if {$proj ne "gnomonic"} {
	return {}
    }

    set radesys [AsdfCelestialRadesys $yamltext]
    if {$radesys eq {}} {
	return {}
    }

    set crpix [AsdfInlineNdarray $body crpix]
    set crval [AsdfInlineNdarray $body crval]
    set cdelt [AsdfInlineNdarray $body cdelt]
    set pc [AsdfInlineNdarray $body pc]
    if {[llength $crpix] != 2 || [llength $crval] != 2 ||
	[llength $cdelt] != 2 || [llength $pc] != 4} {
	return {}
    }

    # Two conversions, and both matter:
    #
    #  - CRPIX is 1-based in FITS but 0-based in gwcs (the node's own
    #    bounding_box runs [-0.5, n-0.5], which is the 0-based range), so
    #    add 1. Off by one here puts the whole image one pixel out.
    #  - FITS CDi_j folds cdelt into the matrix: CD_ij = cdelt_i * pc_ij.
    #    Emitting CDELTi + PCi_j instead would be just as valid FITS, but
    #    the CD form is what every reader handles without argument.
    lassign $crpix cx cy
    lassign $crval rv1 rv2
    lassign $cdelt d1 d2
    lassign $pc p11 p12 p21 p22

    set cards {}
    append cards [AsdfFitsCard CRPIX1 [expr {$cx + 1}]]
    append cards [AsdfFitsCard CRPIX2 [expr {$cy + 1}]]
    append cards [AsdfFitsCard CRVAL1 $rv1]
    append cards [AsdfFitsCard CRVAL2 $rv2]
    append cards [AsdfFitsCard CTYPE1 {'RA---TAN'}]
    append cards [AsdfFitsCard CTYPE2 {'DEC--TAN'}]
    append cards [AsdfFitsCard CUNIT1 {'deg'}]
    append cards [AsdfFitsCard CUNIT2 {'deg'}]
    append cards [AsdfFitsCard CD1_1 [expr {$d1 * $p11}]]
    append cards [AsdfFitsCard CD1_2 [expr {$d1 * $p12}]]
    append cards [AsdfFitsCard CD2_1 [expr {$d2 * $p21}]]
    append cards [AsdfFitsCard CD2_2 [expr {$d2 * $p22}]]
    append cards [AsdfFitsCard RADESYS "'$radesys'"]
    return $cards
}

# Hand synthesized cards to FitsImage::replaceWCSCards. The sentinel first
# line is what tells replaceWCS() these are the file's own WCS rather than a
# user override, so that resetWCS() re-applies them instead of dropping them
# - the same distinction wcsYaml_ already makes for the GWCS path.
proc AsdfAttachWcsCards {cards} {
    global current

    set ch [file tempfile tmpfn]
    fconfigure $ch -translation binary -encoding utf-8
    puts $ch "#ASDF-FITS-WCS"
    puts -nonewline $ch $cards
    close $ch

    try {
	$current(frame) wcs replace 1 $tmpfn
    } finally {
	file delete -force $tmpfn
    }
}

proc AsdfWcsGridShape {entries} {
    set path [AsdfResolvePath $entries data]
    if {$path eq {}} {
	return {}
    }
    set node [AsdfFindNdarrayPath $entries $path]
    if {$node eq {}} {
	return {}
    }
    lassign $node source datatype byteorder shapelist unsupported
    return $shapelist
}

# Two shapes share a pixel grid if their *trailing two* dimensions agree.
# Trailing rather than exact, because a 2-D WCS legitimately describes each
# plane of a rank-3 stack of the same grid: shape is row-major with the
# fastest axis last, so [nplane, ydim, xdim] and [ydim, xdim] are the same
# grid. That still rejects roman/amp33 ([10, 4096, 128]) and the
# border_ref_pix_* set, whose trailing dims are nothing like [4088, 4088].
proc AsdfSameGrid {a b} {
    if {[llength $a] < 2 || [llength $b] < 2} {
	return 0
    }
    return [expr {[lrange $a end-1 end] eq [lrange $b end-1 end]}]
}

proc AsdfAttachWcs {yamltext} {
    global current

    set ch [file tempfile tmpfn]
    fconfigure $ch -translation binary -encoding utf-8
    puts -nonewline $ch $yamltext
    close $ch

    try {
	$current(frame) wcs replace 1 $tmpfn
    } finally {
	file delete -force $tmpfn
    }
}

# --------------------------------------------------------------------
# Array browser (TODO.md Phase 4, design doc SS8 point 5's second option).
#
# pds9 offers both a "<file>:<path>" text entry and an image browser;
# this is the browser half, for the common case of not knowing what a
# non-Roman file contains. Cheap to build: resolving what is in a file
# needs only the YAML tree, never any block bytes.

# Every enumerated array that could actually be displayed, as
# {path label} pairs, largest first so a file's science arrays sort above
# its small WCS coefficient matrices - a real Roman *_cal.asdf enumerates
# 27 ndarrays, of which 12 are 6x6 polynomial coefficients that nobody
# wants to look at as an image. Rank and datatype are filtered by the same
# rules AsdfLoadArray enforces, so nothing offered here can fail on
# selection.
proc AsdfLoadableArrays {entries} {
    set rows {}
    foreach entry $entries {
	lassign $entry path source datatype byteorder shapelist \
	    unsupported maskscalar bitpix

	if {$unsupported ne {}} {
	    continue
	}

	set rank [llength $shapelist]
	if {$rank != 2 && $rank != 3} {
	    continue
	}

	# an empty bitpix is the reader's own answer for "no BITPIX holds
	# this datatype, and no lossless widening does either", so there
	# is no second copy of that table here
	if {$bitpix eq {}} {
	    continue
	}

	set nn 1
	foreach ss $shapelist {
	    set nn [expr {$nn * $ss}]
	}
	lappend rows [list $nn $path \
			  "$path  \[[join $shapelist {, }]\]  $datatype"]
    }

    set result {}
    foreach row [lsort -integer -decreasing -index 0 $rows] {
	lappend result [list [lindex $row 1] [lindex $row 2]]
    }
    return $result
}

# Asks which array to load, and sets `varname` to its path. Returns 1 to
# proceed, 0 to cancel - the same convention MosaicWCSDialog uses, since
# OpenDialog treats them identically. A file with a single loadable array
# selects it without prompting.
proc AsdfPathDialog {fn varname} {
    upvar $varname var

    if {[catch {AsdfArrays $fn} entries]} {
	Error "[msgcat::mc {Unable to load}] $fn: $entries"
	return 0
    }

    set rows [AsdfLoadableArrays $entries]
    switch -- [llength $rows] {
	0 {
	    Error "[msgcat::mc {ASDF: no loadable arrays found}] $fn"
	    return 0
	}
	1 {
	    set var [lindex [lindex $rows 0] 0]
	    return 1
	}
    }

    set rr [AsdfArrayTreeDialog $entries $rows [msgcat::mc {Select Array}]]
    if {$rr == {}} {
	return 0
    }

    set var $rr
    return 1
}

proc AsdfDef {} {
    global asdftree

    # levels of the array picker's tree shown open; deeper ones start closed
    set asdftree(opendepth) 3
    # same, for the header viewer's tree
    set asdftree(headerdepth) 2
}

# Fills treeview `tv' with the file's loadable arrays, as the tree they sit
# in: only the branches leading to a loadable array appear. `entries' is
# AsdfArrays' result (tree order), `rows' AsdfLoadableArrays' (the loadable
# subset, largest first). Item ids are the paths themselves; array items
# carry the `array' tag and their shape and datatype. A node can be both an
# array and a branch - a masked array's `mask' is itself loadable.
# Returns the path of the largest array.
proc AsdfFillArrayTree {tv entries rows} {
    global asdftree

    set loadable {}
    foreach row $rows {
	dict set loadable [lindex $row 0] 1
    }

    foreach entry $entries {
	lassign $entry path source datatype byteorder shapelist
	if {![dict exists $loadable $path]} {
	    continue
	}

	set parent {}
	set prefix {}
	set depth 0
	foreach part [split $path /] {
	    incr depth
	    if {$prefix eq {}} {
		set prefix $part
	    } else {
		append prefix /$part
	    }
	    if {![$tv exists $prefix]} {
		$tv insert $parent end -id $prefix -text $part \
		    -open [expr {$depth < $asdftree(opendepth)}] \
		    -tags branch
	    }
	    set parent $prefix
	}

	$tv item $path -values [list "\[[join $shapelist {, }]\]" $datatype] \
	    -tags array
    }

    return [lindex [lindex $rows 0] 0]
}

# Modal array picker for an ASDF file. Returns the chosen array's path, or
# {} if cancelled.
proc AsdfArrayTreeDialog {entries rows title} {
    global ed

    set w {.asdftree}
    set ed(ok) 0

    DialogCreate $w $title ed(ok)

    set f [ttk::frame $w.ed]
    ttk::scrollbar $f.yscroll -command "$f.tree yview"
    ttk::scrollbar $f.xscroll -command "$f.tree xview" -orient horizontal
    set tv [ttk::treeview $f.tree \
		-yscroll "$f.yscroll set" \
		-xscroll "$f.xscroll set" \
		-selectmode browse \
		-columns {shape datatype} \
		-height 16 \
	       ]
    $tv heading #0 -text [msgcat::mc {Array}] -anchor w
    $tv heading shape -text [msgcat::mc {Shape}] -anchor w
    $tv heading datatype -text [msgcat::mc {Type}] -anchor w
    $tv column #0 -width 300 -stretch 1
    $tv column shape -width 140 -stretch 0
    $tv column datatype -width 80 -stretch 0
    # branches that are not themselves arrays are not choices
    $tv tag configure branch -foreground gray50

    set first [AsdfFillArrayTree $tv $entries $rows]
    $tv see $first
    $tv selection set [list $first]
    $tv focus $first

    grid $f.tree $f.yscroll -sticky news
    grid $f.xscroll -sticky ew
    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 0 -weight 1

    set f [ttk::frame $w.buttons]
    ttk::button $f.ok -text [msgcat::mc {OK}] \
	-command [list AsdfArrayTreeOK $tv]
    ttk::button $f.cancel -text [msgcat::mc {Cancel}] \
	-command {set ed(ok) 0}
    pack $f.ok $f.cancel -side left -expand true -padx 2 -pady 4

    ttk::separator $w.sep -orient horizontal
    pack $w.buttons $w.sep -side bottom -fill x
    pack $w.ed -side top -fill both -expand true

    # double click chooses an array; on a branch it keeps the treeview's
    # own open/close behavior
    bind $tv <Double-1> [list AsdfArrayTreeOK $tv]
    bind $w <Return> [list AsdfArrayTreeOK $tv]

    DialogWait $w ed(ok) $w.buttons.ok

    set rr {}
    if {$ed(ok)} {
	set rr [lindex [$tv selection] 0]
    }

    destroy $w
    return $rr
}

# Accepts the dialog only when the selection is an array.
proc AsdfArrayTreeOK {tv} {
    global ed

    set sel [lindex [$tv selection] 0]
    if {$sel ne {} && [lsearch [$tv item $sel -tags] array] >= 0} {
	set ed(ok) 1
    }
}

# Splits pds9's "<filename>:<path inside the file>" convention into the
# two parts, returning {filename path} with an empty path when there is no
# suffix. Matching pds9 (github.com/asdf-format/pds9) matters here: it is
# the convention Roman users already type, and keeping it means their
# existing muscle memory and scripts carry over to native loading.
#
# The split is decided by whether the part before the last colon actually
# names an existing file, rather than by pattern-matching the text. That
# is what makes it safe on Windows, where "C:/data/x.asdf" has a colon
# that must not be treated as a path separator - "C" is not a file, so it
# is not split. It also degrades sensibly: if nothing splits, the whole
# string is treated as the filename and the caller's "unable to load"
# error quotes it verbatim.
proc AsdfSplitPath {fn} {
    set colon [string last : $fn]
    if {$colon <= 0} {
	return [list $fn {}]
    }

    set base [string range $fn 0 [expr {$colon - 1}]]
    set path [string range $fn [expr {$colon + 1}] end]
    if {$path eq {} || ![file exists $base]} {
	return [list $fn {}]
    }

    return [list $base $path]
}

# --------------------------------------------------------------------
# Masked arrays -> the FITS null conventions.
#
# asdf represents a numpy masked array as an explicit boolean `mask`
# ndarray beside the data, and the core/ndarray schema also allows a
# scalar `mask` naming the sentinel value outright. Both are resolved in
# C++ (asdfApplyMask in fitsy/asdf.C): an integer array keeps its native
# storage and gets a BLANK keyword, a float array gets NaN written into
# the masked pixels. Nothing is promoted to float, which is the whole
# point - FitsData substitutes NAN only at the getValueFloat() boundary
# and minmax skips blank pixels (tksao/frame/fitsdata.C).
#
# Building the header in C++ is also what retired the in-memory FITS file
# this used to synthesize: fitsy's array-spec grammar has no BLANK
# keyword, so the only way to carry one through the array loader was to
# wrap the pixels in a whole FITS file and hand that over instead.
# FitsAsdf makes its own FitsHead and appends the card.

# The mask cases the reader declines, warned about here rather than left
# silent: the array does load, but without its null values, and that is
# worth saying. Everything else - picking the sentinel, the broadcast
# rules, NaN for floats - happens in C++ and needs nothing from here.
proc AsdfMaskWarn {entries path datatype maskscalar} {
    set masknode [AsdfFindNdarrayPath $entries "$path/mask"]
    if {$masknode eq {} && $maskscalar eq {}} {
	return
    }

    # Only the datatypes with a real FITS BITPIX are eligible. That
    # excludes uint16 (it would need FITS's BZERO unsigned-offset
    # convention) and anything that had to be widened, since its BITPIX
    # is no longer the datatype the mask was written against.
    switch -- $datatype {
	uint8 -
	int16 -
	int32 -
	int64 -
	float32 -
	float64 {}
	default {
	    Warning "[msgcat::mc {ASDF: cannot represent null values for datatype, loading without them}] $datatype"
	    return
	}
    }

    if {$masknode ne {}} {
	lassign $masknode msource mdatatype mbyteorder mshape munsupported
	if {$munsupported ne {}} {
	    Warning "[msgcat::mc {ASDF: unreadable mask array, loading without null values}] $path/mask"
	}
    }
}

# One 80-column FITS card. Values are right-justified in columns 11-30,
# per the standard's fixed format.
proc AsdfFitsCard {keyword value} {
    set card [format "%-8s= %20s" $keyword $value]
    return [format "%-80s" [string range $card 0 79]]
}

# The Open/OpenDialog-facing entry point (see ds9/library/open.tcl's Open
# switch), deliberately mirroring LoadFitsFile's {fn layer mode} signature
# so ASDF dispatches like any other format rather than needing its own
# path through the dialog code.
#
# `mode` has no ASDF meaning and is accepted only to keep that signature:
# for FITS it selects load variants (slice, mosaic flavors) that are all
# FITS-container concepts.
#
# `fn` may carry pds9's "<file>:<path>" suffix to name an array other than
# the default `roman/data` - see AsdfSplitPath.
proc LoadAsdfFile {fn layer mode} {
    lassign [AsdfSplitPath $fn] base path
    if {$path eq {}} {
	set path data
    }

    return [AsdfLoadArray $base $path $layer]
}

# --------------------------------------------------------------------
# Scripting surface (TODO.md Phase 3, design doc SS8 point 5) - the XPA /
# SAMP / command-line entry point, parallel to fits.tcl's ProcessFitsCmd.
# Grammar in ds9/parsers/asdfparser.tac + asdflex.fcl (taccle/fickle, the
# in-tree pure-Tcl generators - nothing here goes near bison/flex).

proc ProcessAsdfCmd {varname iname sock fn} {
    upvar $varname var
    upvar $iname i

    global parse
    set parse(sock) $sock
    set parse(fn) $fn

    asdf::YY_FLUSH_BUFFER
    asdf::yy_scan_string [lrange $var $i end]
    asdf::yyparse
    incr i [expr $asdf::yycnt-1]
}

# Unlike FitsCmdLoad there is no socket/stdin variant. The container
# reader resolves an ndarray's `source: N` against the block index written
# at the *end* of the file, so it needs the whole file addressable, not a
# forward-only byte stream - a piped ASDF would have to be spooled to a
# temp file first and gains nothing over naming the file directly. Say so
# rather than silently doing nothing when handed no filename.
proc AsdfCmdLoad {param layer} {
    if {$param == {}} {
	Error [msgcat::mc {ASDF: a filename is required}]
	return
    }

    LoadAsdfFile $param $layer {}
    FinishLoad
}

# Loads one ndarray from an ASDF file into the current frame, through the
# frame's own `load asdf' command - the C++ container reader in
# fitsy/asdf.C, reached the same way `load fits mmap' reaches the FITS
# one. `key` is either a full root-relative path (roman/dq, or any deeper
# path the enumerator reports) or a bare name, resolved by
# AsdfResolvePath - default "data", i.e. Roman's fixed science array.
# Returns 1 on success, 0 on failure (matching ProcessLoad).
proc AsdfLoadArray {fn {key data} {layer {}}} {
    global current
    global loadParam

    if {$current(frame) == {}} {
	CreateFrame
    }

    if {[catch {AsdfArrays $fn} entries]} {
	Error "[msgcat::mc {Unable to load}] $fn: $entries"
	return 0
    }
    if {[catch {AsdfTree $fn} tree]} {
	Error "[msgcat::mc {Unable to load}] $fn: $tree"
	return 0
    }

    set path [AsdfResolvePath $entries $key]
    if {$path == {}} {
	Error "[msgcat::mc {ASDF: ambiguous or unknown array}] $key"
	return 0
    }
    set node [AsdfFindNdarrayPath $entries $path]
    if {$node == {}} {
	Error "[msgcat::mc {ASDF: could not find ndarray}] $path"
	return 0
    }
    lassign $node source datatype byteorder shapelist unsupported \
	maskscalar bitpix

    # Everything a user can get wrong about an array is checked here, not
    # left to the C++ reader, even though it enforces the same rules
    # itself - it has to, since it is what actually reads the block.
    # ProcessLoad discards the frame command's own error text, so a
    # refusal raised down there would reach the user as a bare "Unable to
    # load"; these are the messages anyone actually sees.
    if {[lsearch -exact $unsupported byteorder] >= 0} {
	Error "[msgcat::mc {ASDF: ndarray has no byteorder}] $path"
	return 0
    }
    if {$unsupported ne {}} {
	Error "[msgcat::mc {ASDF: unsupported ndarray view}] $path ([join $unsupported {, }])"
	return 0
    }
    if {$bitpix == {}} {
	Error "[msgcat::mc {ASDF: unsupported ndarray datatype}] $datatype"
	return 0
    }

    # Rank 3 loads as a data cube, which covers roman/amp33
    # ([10,4096,128]) and the border_ref_pix_* set ([10,4096,4] etc).
    # Anything else is refused rather than silently truncated.
    switch -- [llength $shapelist] {
	2 -
	3 {}
	default {
	    Error "[msgcat::mc {ASDF: unsupported ndarray rank}] $path \[[join $shapelist {, }]\]"
	    return 0
	}
    }

    AsdfMaskWarn $entries $path $datatype $maskscalar

    # What the Info panel, the filename display and `xpaget ds9 file`
    # show.
    #
    # The label names the array as well as the file, because that is the
    # question a multi-array format actually raises: a Roman *_cal.asdf
    # has 15 loadable arrays of identical shape and WCS, so "which one am
    # I looking at?" is otherwise unanswerable from the UI.
    #
    # It uses the *tail* of the path, not the whole path, and that is not
    # cosmetic: the default `pds9(infobox,filenametype)` is "root base",
    # which runs the name through FitsImage::root() - and that walks back
    # to the last '/' and keeps only what follows. A full `roman/data`
    # label therefore reads as a directory and the panel shows a bare
    # "data", losing the filename that is the whole point. The tail is
    # also still a form AsdfResolvePath accepts on the way in, so the
    # label stays a usable load spec; the fully-qualified path remains
    # visible in the header viewer for the nested cases where the tail is
    # ambiguous.
    set dispname "[file tail $fn]:[file tail $path]"

    # Provenance for Backup rides along in the same array:
    # ProcessLoadSaveParams persists loadParam per frame, so backup.tcl
    # can reload from the ASDF file itself rather than from a FITS
    # conversion that would keep the pixels but drop both the GWCS - it
    # has no FITS-card representation - and the YAML tree. Absolute,
    # because a backup may be written from any working directory.
    set loadParam(file,type) asdf
    set loadParam(file,mode) {}
    set loadParam(load,type) asdf
    set loadParam(file,name) $dispname
    set loadParam(load,layer) $layer
    set loadParam(asdf,file) [file normalize $fn]
    set loadParam(asdf,path) $path

    set rr [ProcessLoad]

    # Everything below is about the frame's *image*, so none of it applies
    # to a mask: a mask is an overlay, and the frame keeps reporting the
    # image as its file (which is correct). Caching the mask's tree would
    # make the header viewer show the mask product instead of the
    # displayed one - or, if the image came from FITS, replace a real FITS
    # header with an ASDF tree - and attaching the mask's WCS would move
    # the displayed image's coordinates. Both were live bugs: a mask from
    # a same-grid file really did overwrite the frame's WCS.
    set isimage [expr {$layer eq {}}]

    # Cache the YAML tree for the header viewer. Deliberately after
    # ProcessLoad, not before: ProcessLoad clears any previously cached
    # tree for this frame (that is where every format's load funnels
    # through), so setting it first would just be wiped.
    if {$rr && $isimage} {
	AsdfSetTree $current(frame) $fn $tree
    }

    # WCS attachment: best-effort, never fatal - the pixel data is
    # already loaded at this point and stays loaded whatever happens
    # here.
    #
    # Three outcomes, deliberately distinguished rather than all
    # swallowed by one catch. Finding no WCS subtree at all is normal and
    # silent: Roman's *_segm.asdf segmentation products genuinely carry
    # none. A subtree that is found but that extraction or AST then
    # rejects (an unsupported tag or schema version, a GWCS shape
    # yamlchan.c does not cover) is different - the frame is still
    # usable, but say so rather than leave the user wondering why there
    # is no WCS. Warning, not Error: it routes to ds9(msg) for XPA/SAMP
    # callers and a non-modal notice in the GUI, which is what a
    # non-fatal condition should do.
    if {$rr && $isimage} {
	if {[catch {AsdfExtractWcsText $tree $fn} yamltext]} {
	    Warning "[msgcat::mc {ASDF: unable to extract WCS, loading without it}] $yamltext"
	} elseif {$yamltext ne {}} {
	    # Two more outcomes, because the file having a usable WCS does
	    # not mean it describes *this* array. Attaching it regardless
	    # is worse than attaching nothing: the coordinates land inside
	    # the science array's footprint and so look entirely plausible
	    # - roman/amp33 (128x4096 reference pixels) read back the same
	    # sky position as roman/data at the same pixel index. Same-grid
	    # siblings (err, dq, var_poisson, chisq, dumo) must keep the
	    # WCS, which is why this is a grid check rather than "only the
	    # science array gets a WCS".
	    #
	    # The two cases get different severities on purpose:
	    #
	    #   - Grid mismatch is silent. It is the normal, expected
	    #     result of loading a reference-pixel array, the user has
	    #     done nothing wrong, and the shapes are right there in the
	    #     array browser. It also keeps a successful load from
	    #     returning XPA$ERROR, which is what Warning does to an
	    #     xpaset caller.
	    #   - Not being able to identify the WCS's grid at all *is*
	    #     surprising, and it is the case where this guard might be
	    #     dropping a WCS that would have been fine - a non-Roman
	    #     file whose science array is not called "data". Say so.
	    set ref [AsdfWcsGridShape $entries]
	    set cards [AsdfFitsWcsImagingCards $yamltext]
	    if {$ref eq {}} {
		Warning "[msgcat::mc {ASDF: cannot tell which array the WCS describes, loading without it}] $path"
	    } elseif {![AsdfSameGrid $shapelist $ref]} {
		# normal: this array is simply not on the WCS's grid
	    } elseif {$cards ne {} && [catch {AsdfAttachWcsCards $cards} msg]} {
		Warning "[msgcat::mc {ASDF: unable to attach WCS, loading without it}] $msg"
	    } elseif {$cards eq {} && [catch {AsdfAttachWcs $yamltext} msg]} {
		Warning "[msgcat::mc {ASDF: unable to attach WCS, loading without it}] $msg"
	    } elseif {![$current(frame) has wcs wcs]} {
		# `wcs replace` does not report failure: when AstYamlChan
		# cannot build a FrameSet - an unknown tag, a schema
		# version past yamlchan.c's MAKE_TEST ceilings, a GWCS
		# shape it has no path for - it returns cleanly and simply
		# leaves the frame with no WCS. Without this check that is
		# completely silent, which is the one outcome the three
		# branches above were written to avoid. Found by
		# WCS_TEST_PLAN I-5.
		#
		# `has wcs wcs` is the right probe: 1 on a working GWCS
		# frame, 0 when none was built. `has wcs alt` is not - it
		# reads 1 in both cases, because replaceWCSYaml sets
		# wcsAltHeader_ regardless (R2 in the test plan).
		Warning "[msgcat::mc {ASDF: AST could not read this WCS, loading without it}] $path"
	    }
	}
    }

    return $rr
}
