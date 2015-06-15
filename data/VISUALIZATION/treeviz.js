/*
Constructor for the tree drawer. The arguments are provided as a an object, with the
following keys:
    svg - the SVG root (document) element to which all drawn elements are appended
    pptu - pixels per time unit, e.g. if the time units are in MY, how many pixels is 1 MY
    nowx - where the present is on the X axis
    maxWidth - maximum width for any branch, default is 10 pixels
*/
var TreeDrawer = function(args) {
    this.svg = args.svg;
    this.pptu = args.pptu;
    this.nowx = args.nowx;
    this.maxWidth = args.maxWidth || 10;
    this.doc = document;
    this.NS_SVG = 'http://www.w3.org/2000/svg';
    this.NS_XLINK = 'http://www.w3.org/1999/xlink';
    this.NS_XHTML = 'http://www.w3.org/1999/xhtml';
    var NCBI = 'http://ncbi.nlm.nih.gov/';
    this.NCBI_TAXONOMY = NCBI + 'taxonomy/';
    this.NCBI_NUCCORE = NCBI + 'nuccore/';
};

/*
Once the tree drawer is constructed, this is the only method that you would use from the
outside. The argument tree is a JSON object as produced by Bio::Phylo::to_js
*/
TreeDrawer.prototype.drawTree = function(tree) {
    this.markers = tree.getBackboneMarkers();
    this.clades = tree.getCladeMarkers();
    this.maxMarkers = Object.keys(this.markers).length;
    this.recursiveDraw(tree.getRoot());
};

/*
Starting from the root, i.e. invoked as recursiveDraw(root,null), traverses
the tree to draw it, invoking all relevant item drawers.  Age ranges, branches,
and node handlers are invoked in pre-order, the node drawer is invoked in post-order
so that the node lands on top of age ranges and abutting branches.
@param node - the focal node and its attributes to draw
@param parent - optionally, the focal node's parents to draw branches to
 */
TreeDrawer.prototype.recursiveDraw = function(node,parent) {
    var td = this;
    var children = node.getChildren();
    var childCount = children.length;

    // draw the age range for interior nodes
    if ( childCount != 0 ) {
        this.drawAgeRange(node);
    }

    // if not root, draw the branch to the parent
    if ( parent != null ) {
        this.drawBranch(node, parent);
    }

    // if node has a label, draw it. make bold if exemplar. add link.
    if ( node.getName() != null ) {
        var label = this.drawNodeLabel(node);
        var content = this.createNodeContent(node);
        label.onclick = function () {
            td.drawTable(node.getX(),node.getY(),content);
        }
    }

    // recurse
    for ( var i = 0; i < childCount; i++ ) {
        this.recursiveDraw(children[i],node);
    }

    // draw the node if internal
    if ( childCount != 0 ) {
        this.drawNode(node);
    }
};

/*
Draws the age range around a node. This method expects that the node is annotated
with the phylomap annotations map:x and map:y (for the coordinates) and the figtree
annotations fig:height_95_HPD_min and fig:height_95_HPD_max for the age ranges in
some time unit that is convertible to phylomap coordinates by multiplying with the
pixels per time unit (pptu) constant and subtracting that from the x coordinate of
the present (nowx).
 */
TreeDrawer.prototype.drawAgeRange = function(node) {

    // coordinates of the node
    var nx = node.getX();
    var ny = node.getY();

    // have age range
    if ( node.hasAgeRange() ) {

        // rightmost limit x coordinate
        var min_x = this.nowx - ( node.getMinAge() * this.pptu );
        var min_width = Math.abs(min_x - nx);

        // leftmost limit
        var max_x = this.nowx - ( node.getMaxAge() * this.pptu );
        var max_width = Math.abs(max_x - nx);

        var fadeRight = this.createSVGElt('rect',{
            'x'      : nx,
            'y'      : ny - ((this.maxWidth/2)),
            'width'  : min_width,
            'height' : this.maxWidth,
            'style'  : this.makeStyle({ 'fill' : 'url(#fadeRight)', 'stroke' : 'none' })
        });
        this.svg.appendChild(fadeRight);
        var fadeLeft = this.createSVGElt('rect',{
            'x'      : max_x,
            'y'      : ny - ((this.maxWidth/2)),
            'width'  : max_width,
            'height' : this.maxWidth,
            'style'  : this.makeStyle({ 'fill' : 'url(#fadeLeft)', 'stroke' : 'none' })
        });
        this.svg.appendChild(fadeLeft);
    }
};

/*
Draws the focal node, i.e. draws the node glyph (a circle) and any attached behaviours
(i.e. click events that pop up additional node metadata). This method expects that the
node object is annotated with the phylomap predicates map:x, map:y, map:radius and map:branch_width.
Optionally, the node may also have an annotation to indicate that it is the root of a
decomposed clade, using the fig:clade annotation, and it may have the fossil attribute to indicate
that this node was calibrated.
 */
TreeDrawer.prototype.drawNode = function(node) {
    var td = this;
    var nx = node.getX();
    var ny = node.getY();
    var strokeColor = node.getCladeName() ? 'lime' : 'black';
    var nodeColor = node.getFossil() ? 'red' : 'white';
    var circleElt = this.drawCircle(nx, ny, node.getRadius(),{
        'fill'         : nodeColor,
        'stroke'       : strokeColor,
        'stroke-width' : node.getBranchWidth(),
        'cursor'       : 'pointer'
    });
    var content = this.createNodeContent(node);
    circleElt.onclick = function () {
        td.drawTable(nx, ny, content);
    };

    // the element and the node object hold references to each other
    node['node'] = circleElt;
    circleElt.node = node;
};

/*
Draws the branch between the focal node and its parent. The nodes are expected to have
the phylomap attributes map:x, map:y and map:branch_color. The branch is further decorated
to indicate the amount of backbone marker support.
 */
TreeDrawer.prototype.drawBranch = function(node,parent){
    var nx = node.getX();
    var ny = node.getY();
    var px = parent.getX();
    var py = parent.getY();
    var lines = [];

    // calculate branch width (proportional to
    // number of markers
    var width = 2;
    var bbmarkers = node.getBackboneMarkers();
    if (bbmarkers) {
        var mc = Object.keys(bbmarkers).length;
        width = this.maxWidth * ( mc / this.maxMarkers );
    }

    // draw styled lines, one vertical, one horizontal
    var branchStyle = {
        'stroke': this.rgb(node.getBranchColor()),
        'stroke-width': width,
        'stroke-linecap': 'round'
    };
    this.drawLine(px, py, px, ny, branchStyle);
    this.drawLine(px, ny, nx, ny, branchStyle);

    // overlay black dotted line if backbone
    if (bbmarkers) {
        branchStyle['stroke'] = 'black';
        branchStyle['stroke-dasharray'] = '5, ' + ( width + 5 );
        branchStyle['stroke-linecap'] = 'square';
        lines.push(this.drawLine(px, py, px, ny, branchStyle));
        lines.push(this.drawLine(px, ny, nx, ny, branchStyle));
    }

    // the node holds references to the in-pointing branch segments
    node['branch'] = lines;
};

/*
Draws the focal node's label. For the @param node, expects the following
phylomap attributes:
- map:x
- map:y
- map:font_color (optional)
- map:text_horiz_offset (optional)
- map:text_vert_offset (optional)
- map:font_face (optional)
- map:font_size (optional)
- map:font_style (optional)
In addition, optionally expects the attribute 'exemplar'
 */
TreeDrawer.prototype.drawNodeLabel = function(node) {
    var nx = node.getX();
    var ny = node.getY();
    var fontWeight = node.isExemplar() ? 'bold' : 'normal';
    return this.drawText(
        nx + node.getTextHorizOffset(),
        ny + node.getTextVertOffset(),
        node.getName(),
        {
            'stroke': node.getFontColor(),
            'font-family': node.getFontFace(),
            'font-size': node.getFontSize(),
            'font-style': node.getFontStyle(),
            'cursor': 'pointer',
            'font-weight': fontWeight
        }
    );
};

/*
Draws a line glyph.
@param x1 - the starting point x coordinate
@param y1 - the starting point y coordinate
@param x2 - the endpoint x coordinate
@param y2 - the endpoint y coordinate
@param style - a simple object that contains CSS styles
 */
TreeDrawer.prototype.drawLine = function(x1,y1,x2,y2,style) {
    var lineElt = this.createSVGElt('line',{
        'x1'    : x1,
        'y1'    : y1,
        'x2'    : x2,
        'y2'    : y2,
        'style' : this.makeStyle(style)
    });
    this.svg.appendChild(lineElt);
    return lineElt;
};

/*
Draws a circle glyph.
@param cx - the x coordinate of the circle's center
@param cy - the y coordinate of the circle's center
@param r - the radius of the circle
@param style - a simple object that contains CSS styles
 */
TreeDrawer.prototype.drawCircle = function(cx,cy,r,style) {
    var nodeElt = this.createSVGElt('circle',{
        'cx'    : cx,
        'cy'    : cy,
        'r'     : r,
        'style' : this.makeStyle(style)
    });
    this.svg.appendChild(nodeElt);
    return nodeElt;
};

/*
Draws SVG text.
@param x - the x coordinate of the start of the string
@param y - the y coordinate of the start of the string
@param text - the text string to draw
@param style - a simple object that contains CSS styles
@param url - a URL to turn the text into a clickable link
 */
TreeDrawer.prototype.drawText = function(x,y,text,style,url) {
    if ( url ) {
        var aElt = this.doc.createElementNS(this.NS_SVG,'a');
        aElt.setAttributeNS( this.NS_XLINK, 'xlink:href', url );
        aElt.setAttributeNS( this.NS_XLINK, 'xlink:show', 'new');
        this.svg.appendChild(aElt);
    }

    var txtElt = this.doc.createElementNS(this.NS_SVG,'text');
    var txt = this.doc.createTextNode(text);
    txtElt.appendChild(txt);
    txtElt.setAttributeNS( null, 'x', x );
    txtElt.setAttributeNS( null, 'y', y );

    if ( style ) {
        txtElt.setAttributeNS( null, 'style', this.makeStyle(style) );
    }
    if ( url ) {
        aElt.appendChild(txtElt);
        return aElt;
    }
    else {
        this.svg.appendChild(txtElt);
        return txtElt;
    }
};

/*
Draws an (xhtml) table row.
@param row - an object that can contain the following: 'url', a link
             to make clickable. Any values that are functions, which
             are applied as event handlers. Any other, textual, key/
             value pairs are written in the left and right row cells,
             respectively.
@param tr - an xhtml tr (table row) element into which to insert the cells
 */
TreeDrawer.prototype.drawTableRow = function(row,tr) {
    var events = {};
    var textElt;

    // take the URL annotation
    var url = row['url'];
    delete row['url'];

    // take the marker annotation
    var marker = row['marker'];
    delete row['marker'];

    for (var prop in row) {
        if (row.hasOwnProperty(prop)) {

            // item is an event handler
            if (typeof row[prop] === 'function') {
                events[prop] = row[prop];
            }

            // item is textual key/value pair
            else {
                textElt = this.createElt('td',null,this.NS_XHTML);
                textElt.appendChild(this.doc.createTextNode(prop));
                tr.appendChild(textElt);

                var value = this.createElt('td',null,this.NS_XHTML);
                var txt = this.doc.createTextNode(row[prop]);
                if ( url ) {
                    var a = this.createElt('a',{'href':url},this.NS_XHTML);
                    a.appendChild(txt);
                    value.appendChild(a);
                }
                else {
                    value.appendChild(txt);
                }
                tr.appendChild(value);
            }
        }
    }

    // have event handlers, attach them to the itemTitle
    if (Object.keys(events).length > 0) {
        for (var eventType in events) {
            if (events.hasOwnProperty(eventType)) {
                textElt.addEventListener(eventType, events[eventType]);
            }
        }
    }

    // have a URL or events, style itemTitle as a link
    if (Object.keys(events).length > 0) {
        this.addClass(textElt,'link');
    }

    // have a marker ID, attach to element title
    if (marker != null) {
        textElt.setAttributeNS(null, 'title', marker);
    }
};

/*
Draws an xhtml table with node metadata.
@param x - the x coordinate of the table (top-left)
@param y - the y coordinate of the table (top-left)
@param content - an array of node metadata items
 */
TreeDrawer.prototype.drawTable = function(x,y,content) {
    var lines    = content.length;
    var height   = 15;
    var sections = this.countSections(content);

    // probably want to reduce the amount of magic numbers here
    var fo = this.createSVGElt('foreignObject',{
        x      : x,
        y      : y,
        width  : 220,
        height : height * ( lines + 2 ) + sections * ( height + 15 ) + 60
    });
    this.svg.appendChild(fo);

    // create XHTML body, append to foreignObject
    var body = this.createElt('body',null,this.NS_XHTML);
    fo.appendChild(body);

    // create XHTML table, append to body
    var table = this.createElt('table',{
        cellpadding : 0,
        cellspacing : 0,
        style       :'display:none'
    },this.NS_XHTML);
    body.appendChild(table);

    // create XHTML table header, append to table
    var header = this.createElt('th',{ colspan : 2, class : 'widget' },this.NS_XHTML);
    table.appendChild(header);

    // create XHTML close button, insert in header
    var button = this.createElt('button',{ title :'close'},this.NS_XHTML);
    header.appendChild(button);
    button.onclick = function(){ $(fo).fadeOut(400) };

    // append XHTML table rows
    for ( var i = 0; i < lines; i++ ) {
        var tr = this.createElt('tr',null,this.NS_XHTML);

        // line is a single string, use as section header
        if ( typeof content[i] === 'string' ) {
            var th = this.createElt('th',{
                colspan : 2,
                class   : 'sectionHeader'
            },this.NS_XHTML);
            var txt = this.doc.createTextNode(content[i]);
            th.appendChild(txt);
            tr.appendChild(th);
        }

        // line is an item with extra annotations (events? links? key/value data?)
        else if ( content[i] != null ) {
            this.drawTableRow(x, y, content[i], tr);
        }
        table.appendChild(tr);
    }
    $(table).fadeIn(400);
};

/*
XXX use jQuery for this
Adds the provided class name to the provided element.
@param elt - the element to attach the class name to
@param className - the class name to attach to the element
 */
TreeDrawer.prototype.addClass = function(elt,className) {
    var current = elt.getAttributeNS(null, 'class');
    if ( current ) {
        className = current + ' ' + className;
    }
    elt.setAttributeNS(null, 'class', className);
};

/*
XXX use jQuery for this
Removes the provided class name from the provided element
@param elt - the element to remove the class name from
@param className - the class name to remove from the element
 */
TreeDrawer.prototype.removeClass = function(elt,className) {
    var current = elt.getAttributeNS(null, 'class');
    if ( current ) {
        var classes = current.split(' ');
        var spliced = [];
        var numClasses = classes.length;
        for ( var i = 0; i < numClasses; i++ ) {
            if ( classes[i] != className ) {
                spliced.push(classes[i]);
            }
        }
        elt.setAttributeNS(null, 'class', spliced.join(' '));
    }
};

/*
XXX use jQuery for this
Checks if the provided element belongs to the provided class
@param elt - the element to check
@param className - the class name to check
 */
TreeDrawer.prototype.hasClass = function(elt,className) {
    var current = elt.getAttributeNS(null, 'class');
    if ( current ) {
        var classes = current.split(' ');
        var numClasses = classes.length;
        for ( var i = 0; i < numClasses; i++ ) {
            if ( classes[i] == className ) {
                return true;
            }
        }
    }
    return false;
};

/*
XXX maybe make this a private static method
Turns the provided array of three numbers (should be integers between
0 and 255) into a CSS-compatible rgb() statement.
@param triple - an array of three numbers
 */
TreeDrawer.prototype.rgb = function(triplet) {
    return 'rgb(' + triplet.join() + ')';
};

/*
XXX maybe make this a private static method
Turns the provided object with CSS styles into a string.
@param style - a simple object whose keys and values are CSS
syntax (e.g. { color: 'black' } )
 */
TreeDrawer.prototype.makeStyle = function(style){
    var string = '';
    for ( var property in style ) {
        if ( style.hasOwnProperty(property) ) {
            string += property + ':' + style[property] + ';';
        }
    }
    return string;
};

/*
XXX maybe make this a private static method
Shortens the provided string to at most 14 characters, last
three are ellipsis.
@param string - the string to shorten
 */
TreeDrawer.prototype.shorten = function (string) {
    if ( string.length > 11 ) {
        var shortened = string.substring(0,11);
        shortened += '...';
        string = shortened;
    }
    return string;
};

/*
XXX maybe make this a private static method
Counts the number of subsections in a collection of node metadata.
@param content - an array of node metadata
 */
TreeDrawer.prototype.countSections = function(content){
    var count = 0;
    var cl = content.length;
    for ( var i = 0; i < cl; i++ ) {
        if ( typeof content[i] === 'string' ) {
            count++;
        }
    }
    return count;
};

/*
Recursively applies or removes the class name 'painted' to the focal node and
all descendants that are supported by the provided marker. Optionally applies
a callback to each applicable node.
@param node - the focal node to paint
@param marker - if the focal node is supported by this marker (un)paint it
@param remove - boolean to switch between adding or removing the class name 'painted'
@param func - optional callback to apply to the node
 */
TreeDrawer.prototype.recursivePaint = function(node,marker,remove,func) {
    var i = 0;
    var bbmarkers = node.getBackboneMarkers();
    if ( bbmarkers && bbmarkers[marker] ) {
        if ( func ) {
            func(node);
        }

        // if boolean argument set to true, remove the paint
        // from the branches by removing the 'painted' class
        if ( remove ) {

            // tips don't have a node object reference
            if ( node.node ) {
                this.removeClass(node.node,'painted');
            }
            for ( i = 0; i < node.branch.length; i++ ) {
                this.removeClass(node.branch[i],'painted');
            }
        }
        else {
            if ( node.node ) {
                this.addClass(node.node,'painted');
            }
            for ( i = 0; i < node.branch.length; i++ ) {
                this.addClass(node.branch[i],'painted');
            }
        }
    }
    var childCount = node.children.length;
    for ( i = 0; i < childCount; i++ ) {
        this.recursivePaint(node.children[i],marker,remove,func);
    }
};

/*
Animates the branches leading to the focal node (or turns it off when already
animated). The animation consists of shifting the stroke-dashoffset by one,
ten times per second.
@param node - the node whose branches to animate
 */
TreeDrawer.prototype.branchAnimator = function(node) {
    for ( var i = 0; i < node.branch.length; i++ ) {
        var intervalId = node.branch[i].getAttributeNS(null,'intervalId');

        // animation was activated, kill it
        if ( intervalId ) {
            clearInterval(parseInt(intervalId));
            node.branch[i].removeAttributeNS(null,'intervalId');
        }

        // start a new animation
        else {
            intervalId = setInterval((function(elt){
                return function () {
                    var myElt = elt;
                    var oldVal = myElt.style.getPropertyValue('stroke-dashoffset');
                    if ( oldVal ) {
                        oldVal = parseInt(oldVal);
                    }
                    else {
                        oldVal = 0;
                    }
                    var newVal = Number.toString( ( oldVal - 1 ) % 10 );
                    myElt.style.setProperty('stroke-dashoffset',newVal,null);
                };
            }(node.branch[i])),100);
            node.branch[i].setAttributeNS(null,'intervalId',intervalId);
        }
    }
};

/*
XXX - maybe use jQuery for this?
Creates the SVG element with the provided name and attributes
@param name - the tag name of the element
@param attributes - the attributes to attach to the element
 */
TreeDrawer.prototype.createSVGElt = function(name,attributes){
    return this.createElt(name,attributes,this.NS_SVG);
};

/*
Creates the XML element with the provided name, attributes and namespace
@param name - the tag name of the element
@param attributes - the attributes to attach to the element
@param ns - the namespace of the element
 */
TreeDrawer.prototype.createElt = function(name,attributes,ns) {
    var elt = this.doc.createElementNS(ns,name);
    if ( attributes != null ) {
        for ( var property in attributes ) {
            if ( attributes.hasOwnProperty(property) && attributes[property] != null ) {
                elt.setAttributeNS(null, property, attributes[property]);
            }
        }
    }
    return elt;
};

/*
Populates the provided content array with metadata about the markers that apply
to the provided node.
@param markerSet - set of markers that applies to the focal node
@param title - node title (XXX remove me)
@param content - growing array of node metadata content
@param markerLookup - lookup table of marker IDs to marker names
@param node - the focal node
 */
TreeDrawer.prototype.createMarkerContent = function(markerSet,title,content,markerLookup,node) {
    content.push(title);
    var td = this;
    for (var property in markerSet) {
        if (markerSet.hasOwnProperty(property)) {
            var row = { 'marker' : property };

            // create row key: marker name(s)
            var concat = this.shorten(markerLookup[property].join(', '));
            row[concat] = markerSet[property] || property;

            // value is the number of sequences per marker
            if (typeof markerSet[property] === 'number') {
                row['mouseover'] = function () {
                    var marker = this.getAttributeNS(null,'title');
                    for ( var i = 0; i < node.children.length; i++ ) {
                        td.recursivePaint(node.children[i],marker,false);
                    }
                };
                row['mouseout'] = function () {
                    var marker = this.getAttributeNS(null,'title');
                    for ( var i = 0; i < node.children.length; i++ ) {
                        td.recursivePaint(node.children[i],marker,true);
                    }
                };
                row['click'] = function () {
                    if ( td.hasClass(this,'painted') ) {
                        td.removeClass(this,'painted');
                    }
                    else {
                        td.addClass(this,'painted');
                    }
                    var marker = this.getAttributeNS(null,'title');
                    for ( var i = 0; i < node.children.length; i++ ) {
                        td.recursivePaint(node.children[i],marker,false,td.branchAnimator);
                    }
                }
            }
            // value is the accession number
            else {
                row['url'] = this.NCBI_NUCCORE + markerSet[property];
            }
            content.push(row);
        }
    }
};

/*
Given a focal node and a growing array of node metadata content, inserts
a section on fossil calibration points.
@param node - the focal node
@param content - the metadata content to expand
 */
TreeDrawer.prototype.createFossilContent = function(node,content) {
    content.push('Fossil');
    content.push({'name':node.getFossilName()});
    content.push({'minimum age' : node.getFossilMinAge()});
    content.push({'maximum age' : node.getFossilMaxAge()});
};

/*
Populates and returns an array of metadata content for the provided node
@param node - the focal node whose metadata is returned
 */
TreeDrawer.prototype.createNodeContent = function(node) {
    // content for the popup table
    var content = [];

    // add table section for fossils
    if ( node.getFossil() ) {
        this.createFossilContent(node, content);
    }

    // add table section for markers
    if ( node.getBackboneMarkers() ) {
        this.createMarkerContent(node.getBackboneMarkers(), 'Backbone markers', content, this.markers, node);
    }
    if ( node.getCladeMarkers() ) {
        this.createMarkerContent(node.getCladeMarkers(), 'Clade markers', content, this.clades, node);
    }

    // add table section for NCBI taxon ID
    var taxonId = node.getId();
    if (taxonId) {
        content.push('NCBI taxonomy');
        content.push({'taxon id': taxonId, 'url': this.NCBI_TAXONOMY + taxonId });
    }

    // add table section highest posterior density intervals, from TreeAnnotator
    content.push('95% HPD intervals');
    var params = [ 'height', 'length', 'dmv1', 'dmv2', 'dmt' ];
    for ( var i = 0; i < params.length; i++ ) {
        var min = 'fig:' + params[i] + '_95_HPD_min';
        var max = 'fig:' + params[i] + '_95_HPD_max';
        if ( node.hasPredicate(min) && node.hasPredicate(max) ) {
            var key = params[i];
            var val = node.getObject(min,2) + ' .. ' + node.getObject(max,2);
            var data = {};
            data[key] = val;
            content.push(data);
        }
    }

    // add support
    if ( node.hasSupport() ) {
        content.push('Support');
        if ( node.getBootstrap() ) {
            content.push({ 'bootstrap' : node.getBootstrap(2) });
        }
        else if ( node.getPosterior() ) {
            content.push({ 'posterior' : node.getPosterior(2) });
        }
    }

    // add clade name
    var cladeName = node.getCladeName();
    if ( cladeName ) {
        content.push('Backbone decomposition');
        content.push({
            'clade ID' : cladeName,
            'url' : cladeName + '/' + cladeName + '-beast-in.xml'
        });
    }

    return content;
};
