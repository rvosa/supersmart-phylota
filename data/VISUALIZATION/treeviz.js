
/**
* Constructs a TreeDrawer object
* @constructor
* @param {Object} args - a command object with the following keys:
* @param {Object} args.svg - the SVG root (document) element to which all drawn elements are appended
* @param {Number} args.pptu - pixels per time unit, e.g. if the time units are in MY, how many pixels is 1 MY
* @param {Number} args.nowx - where the present is on the X axis
* @param {Number} args.maxWidth - maximum width for any branch, default is 10 pixels
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

/**
* Draws the supplied tree object
* @param {Tree} tree - an instance of the Tree class
*/
TreeDrawer.prototype.drawTree = function(tree) {
    this.markers = tree.getBackboneMarkers();
    this.clades = tree.getCladeMarkers();
    this.maxMarkers = Object.keys(this.markers).length;
    this.makeCladeColors(tree);
    this.recursiveDraw(tree.getRoot(),null);
};

/**
 * Applies distinct color for each clade
 * @param {Tree} tree - an instance of the Tree class
 */
TreeDrawer.prototype.makeCladeColors = function(tree) {

    // collect clade names
    var clades = [];
    tree.getRoot().visitDepthFirst({
        pre : function(node) {
            var name = node.getCladeName();
            if (name)
                clades.push(name);
        }
    });

    // shuffle for more constrasting colors
    function shuffle(o){
        for(var j, x, i = o.length; i; j = Math.floor(Math.random() * i), x = o[--i], o[i] = o[j], o[j] = x);
        return o;
    }
    shuffle(clades);

    // make color lookup
    var lookup = {};
    var steps = clades.length;
    for ( var i = 0; i < steps; i++ ) {
        var R = Math.round( ( 255 / steps ) * i );
        var G = Math.abs( R - 127 );
        var B = 255 - R;
        lookup[clades[i]] = [ R, G, B ];
    }

    // apply colors
    var gray = [ 127, 127, 127 ];
    var color = gray;
    tree.getRoot().visitDepthFirst({
        pre : function(node) {
            node.setBranchColor(color);
            var name = node.getCladeName();
            if (name) {
                color = lookup[name];
                node.setNodeColor(color);
            }
        },
        post : function(node) {
            var name = node.getCladeName();
            if (name)
                color = gray;
        }
    });
};

/**
* Starting from the root, i.e. invoked as recursiveDraw(root,null), traverses
* the tree to draw it, invoking all relevant item drawers.  Age ranges, branches,
* and node handlers are invoked in pre-order, the node drawer is invoked in post-order
* so that the node lands on top of age ranges and abutting branches.
* @param {TreeNode} node - the focal node and its attributes to draw
* @param {TreeNode} parent - optionally, the focal node's parents to draw branches to
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

/**
* Draws the age range around a node.
* @param {TreeNode} node - focal node
*/
TreeDrawer.prototype.drawAgeRange = function(node) {

    // coordinates of the node
    var nx = node.getX();
    var ny = node.getY();

    // have age range
    if ( node.hasAgeRange() ) {
        var y = ny - ( this.maxWidth / 2 );

        // rightmost limit x coordinate
        var min_x = this.nowx - ( node.getMinAge() * this.pptu );
        var min_width = Math.abs(min_x - nx);

        // leftmost limit
        var max_x = this.nowx - ( node.getMaxAge() * this.pptu );
        var max_width = Math.abs(max_x - nx);

        var fadeRight = this.createSVGElt('rect',{
            x      : nx,
            y      : y,
            width  : min_width,
            height : this.maxWidth,
            style  : this.makeStyle({ fill : 'url(#fadeRight)', stroke : 'none' })
        });
        this.svg.appendChild(fadeRight);
        var fadeLeft = this.createSVGElt('rect',{
            x      : ( nx - max_width ),
            y      : y,
            width  : max_width,
            height : this.maxWidth,
            style  : this.makeStyle({ fill : 'url(#fadeLeft)', stroke : 'none' })
        });
        this.svg.appendChild(fadeLeft);
    }
};

/**
* Draws the focal node, i.e. draws the node glyph (a circle) and any attached behaviours
* (i.e. click events that pop up additional node metadata).
* @param {TreeNode} node - the focal node
*/
TreeDrawer.prototype.drawNode = function(node) {
    var td = this;
    var nx = node.getX();
    var ny = node.getY();
    var nodeColor = node.getCladeName() ? node.getNodeColor(true) : 'white';
    var strokeColor = node.getFossil() ? 'red' : 'black';
    var circleElt = this.drawCircle(nx, ny, node.getRadius(),{
        fill           : nodeColor,
        stroke         : strokeColor,
        'stroke-width' : node.getBranchWidth(),
        cursor         : 'pointer'
    });
    var content = this.createNodeContent(node);
    circleElt.onclick = function () {
        td.drawTable(nx, ny, content);
    };

    // add support value on node
    if ( node.hasSupport() ) {
        var value;
        if ( node.getBootstrap() ) {
            value = Math.round( node.getBootstrap() * 100 );
        }
        else {
            value = node.getPosterior(2);
        }
        this.drawText(nx - 28,ny - 8,value.toString(),null,null);
    }

    // the element and the node object hold references to each other
    node['node'] = circleElt;
    circleElt.node = node;
};

/**
* Draws the branch between the focal node and its parent.
* @param {TreeNode} node - the focal node
* @param {TreeNode} parent - the focal node's parent
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
        stroke           : this.rgb(node.getBranchColor()),
        'stroke-width'   : width,
        'stroke-linecap' : 'round'
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

/**
* Draws the focal node's label.
* @param {TreeNode} node - the focal node
* @return {SVGElement}
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
            cursor        : 'pointer',
            'font-family' : node.getFontFace(),
            'font-size'   : node.getFontSize(),
            'font-style'  : node.getFontStyle(),
            'font-weight' : fontWeight
        },
        null
    );
};

/**
* Draws a line glyph.
* @param {Number} x1 - the starting point x coordinate
* @param {Number} y1 - the starting point y coordinate
* @param {Number} x2 - the endpoint x coordinate
* @param {Number} y2 - the endpoint y coordinate
* @param {Object} style - a simple object that contains CSS styles
* @return {SVGElement}
*/
TreeDrawer.prototype.drawLine = function(x1,y1,x2,y2,style) {
    var lineElt = this.createSVGElt('line',{
        x1    : x1,
        y1    : y1,
        x2    : x2,
        y2    : y2,
        style : this.makeStyle(style)
    });
    this.svg.appendChild(lineElt);
    return lineElt;
};

/**
* Draws a circle glyph.
* @param {Number} cx - the x coordinate of the circle's center
* @param {Number} cy - the y coordinate of the circle's center
* @param {Number} r - the radius of the circle
* @param {Object} style - a simple object that contains CSS styles
* @return {SVGElement}
*/
TreeDrawer.prototype.drawCircle = function(cx,cy,r,style) {
    var nodeElt = this.createSVGElt('circle',{
        cx    : cx,
        cy    : cy,
        r     : r,
        style : this.makeStyle(style)
    });
    this.svg.appendChild(nodeElt);
    return nodeElt;
};

/**
* Draws SVG text.
* @param {Number} x - the x coordinate of the start of the string
* @param {Number} y - the y coordinate of the start of the string
* @param {String} text - the text string to draw
* @param {Object} style - a simple object that contains CSS styles
* @param {String} url - a URL to turn the text into a clickable link
* @return {SVGElement}
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

/**
* Draws an (xhtml) table row.
* @param {Object} row - an object that can contain the following: 
* @param {String} row.url - a link to make clickable. 
*   Any values that are functions, which are applied as event handlers. 
*   Any other, textual, key/value pairs are written in the left and right row cells, respectively.
* @param {Element} tr - an xhtml tr (table row) element into which to insert the cells
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

/**
* Draws an xhtml table with node metadata.
* @param {Number} x - the x coordinate of the table (top-left)
* @param {Number} y - the y coordinate of the table (top-left)
* @param {Object[]} content - an array of node metadata items
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
    var args = {
        title : 'close',
        type  : 'image',
        src   : "data:image/gif;base64,R0lGODlhEwATAOYAAAAAAP////f399zc3NXV1dDQ0MzMzMXFxcPDw8LCwr+/v76+vru7u7m5uba2trS0tLOzs7KysrGxsbCwsK+vr66urq2traysrKurq6qqqqmpqaioqKenp6ampqWlpaSkpKOjo6KioqGhoaCgoJ+fn56enp2dnZycnJubm5qampmZmZiYmJeXl5aWlpWVlZSUlJOTk5KSkpGRkZCQkI+Pj46Ojo2NjYyMjIuLi4qKiomJiYiIiIeHh4aGhoWFhYSEhIODg4KCgoGBgX9/f35+fn19fXx8fHt7e3p6enl5eXh4eHd3d3Z2dnV1dXR0dHNzc3FxcW9vb25ubmxsbGlpaWhoaGVlZWRkZGJiYmFhYWBgYF9fX15eXl1dXVtbW1paWllZWVhYWFdXV1VVVf///wAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAGQALAAAAAATABMAAAf/gAoUHyUqLTAvMjI0NDk5Oj9FQUkSCgkNExcaHB4gICMlKCwyOklPHAgHDBAVGBsdHh4hIycrMDdFWyMMCw4VGR2fJigrLIk1OURYJxEPFQMBIyYqLjE1NgEFOjxDVioXFhkB49TWOOMBPD3dLhsZHCDoMDM36D89PURXLh8eHyIm0NVAByQIkB/6aIwIMeuECnTjfAARYrDIFRsnSIwg4bAFuh0SKQo5giXHChMlSjiEGKDHDyFChiTJ0uNFihPD0MUYJ+DlECJMtASZ0SIFCnSJaKCLecRJFyI3YLBA5+LYuXFFkkDxcoRHDRgEAryAscjGtgAGjiyRAkZJEB01RmbMYFTjhg6XQrI6mRLmyZEgPXTk2LFDnU8jSZpAqTIGSZQnSpAgSUJZiZImTqBEmUIlixgkbLmI7rJFy5Ytor18+SKmdSAAOw=="
    };
    var button = this.createElt('input',args,this.NS_XHTML);
    header.appendChild(button);
    button.onclick = function(){ $(fo).fadeOut(400) };

    // append XHTML table rows
    for ( var i = 0; i < lines; i++ ) {
        var tr = this.createElt('tr',null,this.NS_XHTML);

        // line is a single string, use as section header
        if ( typeof content[i] === 'string' ) {
            var th = this.createElt('th',{
                colspan : 2,
                'class' : 'sectionHeader'
            },this.NS_XHTML);
            var txt = this.doc.createTextNode(content[i].toString());
            th.appendChild(txt);
            tr.appendChild(th);
        }

        // line is an item with extra annotations (events? links? key/value data?)
        else if ( content[i] != null ) {
            this.drawTableRow(content[i], tr);
        }
        table.appendChild(tr);
    }
    $(table).fadeIn(400);
};

/**
* XXX use jQuery for this
* Adds the provided class name to the provided element.
* @param {Element} elt - the element to attach the class name to
* @param {String} className - the class name to attach to the element
*/
TreeDrawer.prototype.addClass = function(elt,className) {
    var current = elt.getAttributeNS(null, 'class');
    if ( current ) {
        className = current + ' ' + className;
    }
    elt.setAttributeNS(null, 'class', className);
};

/**
* XXX use jQuery for this
* Removes the provided class name from the provided element
* @param {Element} elt - the element to remove the class name from
* @param {String} className - the class name to remove from the element
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

/**
* XXX use jQuery for this
* Checks if the provided element belongs to the provided class
* @param {Element} elt - the element to check
* @param {String} className - the class name to check
* @return {Boolean}
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

/**
* XXX maybe make this a private static method
* Turns the provided array of three numbers (should be integers between
* 0 and 255) into a CSS-compatible rgb() statement.
* @param {Number[]} triplet - an array of three numbers
* @return {String}
 */
TreeDrawer.prototype.rgb = function(triplet) {
    return 'rgb(' + triplet.join() + ')';
};

/**
* XXX maybe make this a private static method
* Turns the provided object with CSS styles into a string.
* @param {Object} style - a simple object whose keys and values are CSS
* @return {String}
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

/**
* XXX maybe make this a private static method
* Shortens the provided string to at most 14 characters, last
* three are ellipsis.
* @param {String} string - the string to shorten
* @return {String}
*/
TreeDrawer.prototype.shorten = function (string) {
    if ( string.length > 11 ) {
        var shortened = string.substring(0,11);
        shortened += '...';
        string = shortened;
    }
    return string;
};

/**
* XXX maybe make this a private static method
* Counts the number of subsections in a collection of node metadata.
* @param {Object[]} content - an array of node metadata
* @return {Number}
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

/**
* Recursively applies or removes the class name 'painted' to the focal node and
* all descendants that are supported by the provided marker. Optionally applies
* a callback to each applicable node.
* @param {TreeNode} node - the focal node to paint
* @param {String} marker - if the focal node is supported by this marker (un)paint it
* @param {Boolean} remove - boolean to switch between adding or removing the class name 'painted'
* @param {Function} func - optional callback to apply to the node
*/
TreeDrawer.prototype.recursivePaint = function(node,marker,remove,func) {
    var td = this;
    node.visitDepthFirst({
        pre : function(n) {
            var backboneMarkers = n.getBackboneMarkers();
            if ( backboneMarkers && backboneMarkers[marker] ) {

                // optionally used to activate/end the branch animator
                if ( func ) {
                    func(n);
                }

                // if boolean argument set to true, remove the paint
                // from the branches by removing the 'painted' class
                if ( remove ) {

                    // tips don't have a node object reference
                    if ( n['node'] ) {
                        td.removeClass(n['node'],'painted');
                    }
                    for ( i = 0; i < n['branch'].length; i++ ) {
                        td.removeClass(n['branch'][i],'painted');
                    }
                }
                else {
                    if ( n['node'] ) {
                        td.addClass(n['node'],'painted');
                    }
                    for ( i = 0; i < n['branch'].length; i++ ) {
                        td.addClass(n['branch'][i],'painted');
                    }
                }
            }

        }
    });
};

/**
* Animates the branches leading to the focal node (or turns it off when already
* animated). The animation consists of shifting the stroke-dashoffset by one,
* ten times per second.
* @param {TreeNode} node - the node whose branches to animate
*/
TreeDrawer.prototype.branchAnimator = function(node) {
    for ( var i = 0; i < node['branch'].length; i++ ) {
        var intervalId = node['branch'][i].getAttributeNS(null,'intervalId');

        // animation was activated, kill it
        if ( intervalId ) {
            clearInterval(parseInt(intervalId));
            node['branch'][i].removeAttributeNS(null,'intervalId');
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
                    var dashArray = myElt.style.getPropertyValue('stroke-dasharray').split(' ');
                    var newVal = ( oldVal - 1 ) % dashArray[1];
                    var newValString = newVal.toString();
                    myElt.style.setProperty('stroke-dashoffset',newValString,null);
                };
            }(node['branch'][i])),100);
            node['branch'][i].setAttributeNS(null,'intervalId',intervalId);
        }
    }
};

/**
* XXX - maybe use jQuery for this?
* Creates the SVG element with the provided name and attributes
* @param {String} name - the tag name of the element
* @param {Object} attributes - the attributes to attach to the element
* @return {Element}
*/
TreeDrawer.prototype.createSVGElt = function(name,attributes){
    return this.createElt(name,attributes,this.NS_SVG);
};

/**
* Creates the XML element with the provided name, attributes and namespace
* @param {String} name - the tag name of the element
* @param {Object} attributes - the attributes to attach to the element
* @param {String} ns - the namespace of the element
* @return {Element}
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

/**
* Populates the provided content array with metadata about the markers that apply
* to the provided node.
* @param {Object} markerSet - set of markers that applies to the focal node
* @param {String} title - node title (XXX remove me)
* @param {Object[]} content - growing array of node metadata content
* @param {Object} markerLookup - lookup table of marker IDs to marker names
* @param {Object} node - the focal node
*/
TreeDrawer.prototype.createMarkerContent = function(markerSet,title,content,markerLookup,node) {
    content.push(title);
    var td = this;
    for (var property in markerSet) {
        if (markerSet.hasOwnProperty(property)) {
            var row = { marker : property };

            // create row key: marker name(s)
            var concat = this.shorten(markerLookup[property].join(', '));
            row[concat] = markerSet[property] || property;

            // value is the number of sequences per marker
            if (typeof markerSet[property] === 'number') {
                row['mouseover'] = function () {
                    var children = node.getChildren();
                    var marker = this.getAttributeNS(null,'title');
                    for ( var i = 0; i < children.length; i++ ) {
                        td.recursivePaint(children[i],marker,false,null);
                    }
                };
                row['mouseout'] = function () {
                    var children = node.getChildren();
                    var marker = this.getAttributeNS(null,'title');
                    for ( var i = 0; i < children.length; i++ ) {
                        td.recursivePaint(children[i],marker,true,null);
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
                    var children = node.getChildren();
                    for ( var i = 0; i < children.length; i++ ) {
                        td.recursivePaint(children[i],marker,false,td.branchAnimator);
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

/**
* Given a focal node and a growing array of node metadata content, inserts
* a section on fossil calibration points.
* @param {Object} node - the focal node
* @param {Object[]} content - the metadata content to expand
*/
TreeDrawer.prototype.createFossilContent = function(node,content) {
    content.push('Fossil');
    content.push({ 'name'        : node.getFossilName()});
    content.push({ 'minimum age' : node.getFossilMinAge()});
    content.push({ 'maximum age' : node.getFossilMaxAge()});
};

/**
* Populates and returns an array of metadata content for the provided node
* @param node - the focal node whose metadata is returned
* @return {Object[]}
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
            content.push({ bootstrap : node.getBootstrap(2) });
        }
        else if ( node.getPosterior() ) {
            content.push({ posterior : node.getPosterior(2) });
        }
    }

    // add clade name
    var cladeName = node.getCladeName();
    if ( cladeName ) {
        content.push('Backbone decomposition');
        content.push({
            'clade ID' : cladeName,
            'url'      : cladeName + '/' + cladeName + '-beast-in.xml'
        });
    }

    return content;
};
