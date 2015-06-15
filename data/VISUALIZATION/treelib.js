
var Tree = function(json) {
    this.root = new TreeNode(json);

    // inner function to wrap the JSON data as objects
    function wrapRecurse(node) {
        var children = node.getChildren();
        var wrapped = [];
        for ( var i = 0; i < children.length; i++ ) {
            var wrappedNode = new TreeNode(children[i],parent);
            wrapped.push(wrappedNode);
            wrapRecurse(wrappedNode);
        }
        node.children = wrapped;
    }
    wrapRecurse(this.root);
};

Tree.prototype.getRoot = function () {
    return this.root;
};

Tree.prototype.getBackboneMarkers = function () {
    return this.root['backbone_markers'];
};

Tree.prototype.getCladeMarkers = function () {
    return this.root['clade_markers'];
};

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

TreeNode.prototype.getX = function() {
    return this['map:x'];
};

TreeNode.prototype.getY = function() {
    return this['map:y'];
};

TreeNode.prototype.getChildren = function () {
    return this['children'];
};

TreeNode.prototype.getParent = function () {
    return this.parent;
};

TreeNode.prototype.getName = function () {
    return this['binomial'];
};

TreeNode.prototype.getMinAge = function () {
    return this['fig:height_95_HPD_min'];
};

TreeNode.prototype.getMaxAge = function () {
    return this['fig:height_95_HPD_max'];
};

TreeNode.prototype.getRadius = function () {
    return this['map:radius'] || 5;
};

TreeNode.prototype.getBranchWidth = function () {
    return this['map:branch_width'] || 2;
};

TreeNode.prototype.getBranchColor = function () {
    return this['map:branch_color'] || 'black';
};

TreeNode.prototype.getTextHorizOffset = function () {
    return this['map:text_horiz_offset'] || 12;
};

TreeNode.prototype.getTextVertOffset = function () {
    return this['map:text_vert_offset'] || 3;
};

TreeNode.prototype.getFontColor = function () {
    return this['map:font_color'] || 'black';
};

TreeNode.prototype.getFontFace = function () {
    return this['map:font_face'] || 'Verdana';
};

TreeNode.prototype.getFontSize = function () {
    return this['map:font_size'] || '10px';
};

TreeNode.prototype.getFontStyle = function () {
    return this['map:font_style'] || 'italic';
};

TreeNode.prototype.getBackboneMarkers = function () {
    return this['backbone_markers'];
};

TreeNode.prototype.getCladeMarkers = function () {
    return this['clade_markers'];
};

TreeNode.prototype.getId = function () {
    return this['guid'];
};

TreeNode.prototype.getBootstrap = function (precision) {
    if ( precision ) {
        parseFloat(this['fig:bootstrap']).toFixed(precision);
    }
    else {
        return this['fig:bootstrap'];
    }
};

TreeNode.prototype.getPosterior = function (precision) {
    if ( precision ) {
        parseFloat(this['fig:posterior']).toFixed(precision);
    }
    else {
        return this['fig:posterior'];
    }
};

TreeNode.prototype.getCladeName = function () {
    return this['fig:clade'];
};

TreeNode.prototype.getObject = function(predicate,precision) {
    if (precision) {
        return parseFloat(this[predicate]).toFixed(precision);
    }
    else {
        return this[predicate];
    }
};

TreeNode.prototype.getFossilName = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['name'];
    }
    return null;
};

TreeNode.prototype.getFossilMinAge = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['min_age'];
    }
    return null;
};

TreeNode.prototype.getFossilMaxAge = function () {
    var fossil = this.getFossil();
    if ( fossil ) {
        return fossil['max_age'];
    }
    return null;
};

TreeNode.prototype.getFossil = function () {
    return this['fossil'];
};

TreeNode.prototype.hasPredicate = function(predicate) {
    return this.hasOwnProperty(predicate)
};

TreeNode.prototype.hasSupport = function () {
    return this.getBootstrap() || this.getPosterior();
};

TreeNode.prototype.hasAgeRange = function () {
    return this.getMinAge() && this.getMaxAge();
};

TreeNode.prototype.isLeaf = function () {
    return this['children'].length == 0;
};

TreeNode.prototype.isRoot = function () {
    return this.parent == null;
};

TreeNode.prototype.isExemplar = function () {
    return this['exemplar'];
};