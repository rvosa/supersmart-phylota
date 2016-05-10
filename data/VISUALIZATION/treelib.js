
/**
 * Constructs a tree object.
 * @param {Object} json - as produced by Bio::Phylo::to_js
 * @constructor
 */

var Tree = function(json) {
    this.root = new TreeNode(json);

    // inner function to wrap the JSON data as objects
    function wrapRecurse(node) {
        var children = node.getChildren();
        var wrapped = [];
        for ( var i = 0; i < children.length; i++ ) {
            var wrappedNode = new TreeNode(children[i],node);
            wrapped.push(wrappedNode);
            wrapRecurse(wrappedNode);
        }
        node.children = wrapped;
    }
    wrapRecurse(this.root);
};

/**
 * Returns the root node of the tree
 * @returns {TreeNode|*}
 */
Tree.prototype.getRoot = function () {
    return this.root;
};

/**
 * Returns a lookup of all the backbone markers in the tree
 * @returns {*}
 */
Tree.prototype.getBackboneMarkers = function () {
    return this.root['backbone_markers'];
};

/**
 * Returns a lookup of all the clade markers in the tree
 * @returns {*}
 */
Tree.prototype.getCladeMarkers = function () {
    return this.root['clade_markers'];
};

/**
 * Constructs a TreeNode object
 * @param {Object} json - as produced by Bio::Phylo::to_js
 * @param {TreeNode} parent - optional parent node
 * @constructor
 */
var TreeNode = function(json,parent) {

    // extend the object with the JSON data properties
    for ( var property in json ) {
        if ( json.hasOwnProperty(property) ) {
            this[property] = json[property];
        }
    }

    // the JSON node object doesn't know about its parent,
    // but this is useful to have so it could be passed
    // separately to the constructor
    this.parent = parent;
};

/**
 * Gets the X-coordinate of the node
 * @returns {Number}
 */
TreeNode.prototype.getX = function() {
    return this['map:x'];
};

/**
 * Gets the Y-coordinate of the node
 * @returns {Number}
 */
TreeNode.prototype.getY = function() {
    return this['map:y'];
};

/**
 * Returns an array of child nodes
 * @returns {TreeNode[]}
 */
TreeNode.prototype.getChildren = function () {
    return this['children'];
};

/**
 * Returns the parent node
 * @returns {TreeNode|*}
 */
TreeNode.prototype.getParent = function () {
    return this.parent;
};

/**
 * Returns a human readable name
 * @returns {String}
 */
TreeNode.prototype.getName = function () {
    return this['binomial'];
};

/**
 * Returns the minimum age of the node
 * @returns {Number}
 */
TreeNode.prototype.getMinAge = function () {
    return this['fig:height_95_HPD_min'];
};

/**
 * Returns the maximum age of the node
 * @returns {Number}
 */
TreeNode.prototype.getMaxAge = function () {
    return this['fig:height_95_HPD_max'];
};

/**
 * Returns the node radius
 * @returns {Number}
 */
TreeNode.prototype.getRadius = function () {
    return this['map:radius'] || 5;
};

/**
 * Returns the branch width
 * @returns {Number}
 */
TreeNode.prototype.getBranchWidth = function () {
    return this['map:branch_width'] || 2;
};

/**
 * Creates an RGB string
 * @param {Number[]} triplet
 * @returns {String}
 */
function makeRGB(triplet) {
    return 'rgb(' + triplet.join() + ')';
}

/**
 * Returns the branch color
 * @param {Boolean} rgb - whether to return an rgb string
 * @returns {String|Number[]}
 */
TreeNode.prototype.getBranchColor = function (rgb ) {
    if ( rgb ) {
        if ( this['map:branch_color'] ) {
            return makeRGB(this['map:branch_color']);
        }
    }
    return this['map:branch_color'] || 'black';
};

/**
 * Returns the node color
 * @param {Boolean} rgb - whether to return an rgb string
 * @returns {String|Number[]}
 */
TreeNode.prototype.getNodeColor = function (rgb ) {
    if ( rgb ) {
        if ( this['map:node_color'] ) {
            return makeRGB(this['map:node_color']);
        }
    }
    return this['map:node_color'] || 'black';
};

/**
 * Returns text horizontal offset
 * @returns {Number}
 */
TreeNode.prototype.getTextHorizOffset = function () {
    return this['map:text_horiz_offset'] || 12;
};

/**
 * Returns text vertical offset
 * @returns {Number}
 */
TreeNode.prototype.getTextVertOffset = function () {
    return this['map:text_vert_offset'] || 3;
};

/**
 * Returns font color
 * @param {Boolean} rgb - whether to return an rgb string
 * @returns {String|Number[]}
 */
TreeNode.prototype.getFontColor = function (rgb) {
    if ( rgb ) {
        if ( this['map:font_color'] ) {
            return makeRGB(this['map:font_color']);
        }
    }
    return this['map:font_color'] || 'black';
};

/**
 * Returns font face
 * @returns {String}
 */
TreeNode.prototype.getFontFace = function () {
    return this['map:font_face'] || 'Verdana';
};

/**
 * Returns font size (px)
 * @returns {String}
 */
TreeNode.prototype.getFontSize = function () {
    return this['map:font_size'] || '10px';
};

/**
 * Returns font style
 * @returns {String}
 */
TreeNode.prototype.getFontStyle = function () {
    return this['map:font_style'] || 'italic';
};

/**
 * Returns list of backbone markers for the focal node
 * @returns {String[]}
 */
TreeNode.prototype.getBackboneMarkers = function () {
    return this['backbone_markers'];
};

/**
 * Returns list of clade markers for the focal node
 * @returns {String[]}
 */
TreeNode.prototype.getCladeMarkers = function () {
    return this['clade_markers'];
};

/**
 * Returns a machine-readable identified
 * @returns {String}
 */
TreeNode.prototype.getId = function () {
    return this['guid'];
};

/**
 * Returns the bootstrap support value
 * @param {Number} precision - number of decimal places
 * @returns {Number}
 */
TreeNode.prototype.getBootstrap = function (precision) {
    if ( precision ) {
        return parseFloat(this['fig:bootstrap']).toFixed(precision);
    }
    else {
        return this['fig:bootstrap'];
    }
};

/**
 * Returns the posterior probability
 * @param {Number} precision - number of decimal places
 * @returns {Number}
 */
TreeNode.prototype.getPosterior = function (precision) {
    if ( precision ) {
        return parseFloat(this['fig:posterior']).toFixed(precision);
    }
    else {
        return this['fig:posterior'];
    }
};

/**
 * Returns the clade identifier
 * @returns {String}
 */
TreeNode.prototype.getCladeName = function () {
    return this['fig:clade'];
};

/**
 * Given any predicate, returns the associated object value
 * @param {String} predicate
 * @param {Number} precision - number of decimal places (optional)
 * @returns {Object}
 */
TreeNode.prototype.getObject = function(predicate,precision) {
    if (precision) {
        return parseFloat(this[predicate]).toFixed(precision);
    }
    else {
        return this[predicate];
    }
};

/**
 * Returns the name of the fossil attached to the node
 * @returns {String}
 */
TreeNode.prototype.getFossilName = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['name'];
    }
    return null;
};

/**
 * Returns the minimum age of the fossil attached to the node
 * @returns {Number}
 */
TreeNode.prototype.getFossilMinAge = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['min_age'];
    }
    return null;
};

/**
 * Returns the maximum age of the fossil attached to the node
 * @returns {Number}
 */
TreeNode.prototype.getFossilMaxAge = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['max_age'];
    }
    return null;
};

/**
 * Returns the fossil object attached to the node
 * @returns {Object}
 */
TreeNode.prototype.getFossil = function () {
    return this['fossil'];
};

/**
 * Returns whether the node is annotated with the predicate
 * @param {String} predicate
 * @returns {Boolean}
 */
TreeNode.prototype.hasPredicate = function(predicate) {
    return this.hasOwnProperty(predicate)
};

/**
 * Returns any support value
 * @returns {Number}
 */
TreeNode.prototype.hasSupport = function () {
    return this.getBootstrap() || this.getPosterior();
};

/**
 * Returns whether the node has an age range
 * @returns {Boolean}
 */
TreeNode.prototype.hasAgeRange = function () {
    return this.getMinAge() && this.getMaxAge();
};

/**
 * Returns whether the node is terminal
 * @returns {Boolean}
 */
TreeNode.prototype.isLeaf = function () {
    return this['children'].length == 0;
};

/**
 * Returns whether the node is the root
 * @returns {Boolean}
 */
TreeNode.prototype.isRoot = function () {
    return this.parent == null;
};

/**
 * Returns whether the node is an exemplar taxon
 * @returns {Boolean}
 */
TreeNode.prototype.isExemplar = function () {
    return this['exemplar'];
};

/**
 * Sets the branch color to the supplied value
 * @param {Number[]} color - an RGB color triplet
 */
TreeNode.prototype.setBranchColor = function (color) {
    this['map:branch_color'] = color;
};

/**
 * Sets the node color to the supplied value
 * @param {Number[]} color - an RGB color triplet
 */
TreeNode.prototype.setNodeColor = function (color) {
    this['map:node_color'] = color;
};

/**
 * Traverses from the focal node in depth-first traversal
 * @param {Object} args - visitor functions to apply
 * @param {Function} args.pre - pre-order function
 * @param {Function} args.post - post-order function
 */
TreeNode.prototype.visitDepthFirst = function(args) {
    if ( args['pre'] ) {
        args['pre'](this);
    }
    var children = this.getChildren();
    for ( var i = 0; i < children.length; i++ ) {
        children[i].visitDepthFirst(args);
    }
    if ( args['post'] ) {
        args['post'](this);
    }
};
