classdef XMLElement < handle
    %DOMNode Convenient representation of an XML DOM node.
    %    This class imitates the key bits of Python's lxml.etree.Element.
    % Example
    %   % parse .xml into an XMLElement object
    %   p = XMLElement.parse('protocol.hpdd');
    %
    %   % find first child element with matching tag
    %   fluids = p.find('Fluids');
    %
    %   % find *all* matching elements at any depth
    %   f_elts = fluids.iter('Fluid');
    %
    %   % print the text content of all fluids' "Name" element
    %   for i = 1:length(f_elts)
    %       disp(f_elts(i).find('Name').text);
    %   end


    properties (SetAccess = protected)
        tag = ''
        attributes = struct
        text = ''
        tail = ''
        children = []
    end

    methods

        function value = get(self, name)
            for i = 1:length(self.attributes)
                if strcmp(self.attributes(i).name, name)
                    value = self.attributes(i).value;
                    return;
                end
            end
            % FIXME is this reasonable, or should we throw?
            value = '';
        end

        function elt = find(self, tag)
            % Note: Only supports simple tag names, not XPath!
            elt = XMLElement.empty;
            for i = 1:length(self.children)
                if strcmp(self.children(i).tag, tag)
                    elt = self.children(i);
                    return
                end
            end
        end

        function elts = iter(self, tag)
            % Note: Returns an array, not an iterator as in lxml.etree.
            if nargin == 1
                tag = '';
            end
            elts = XMLElement.empty;
            if isempty(tag) || strcmp(self.tag, tag)
                elts(end+1) = self;
            end
            for i = 1:length(self.children)
                elts = [elts self.children(i).iter(tag)]; %#ok<AGROW>
            end
        end

    end

    methods (Static)

        function elt = parse(filename)
            elt = parseXML(filename);
        end

    end

end


% ======= Following code based on example in 'xmlread' doc page =======

% ----- Local function PARSEXML -----
function elt = parseXML(filename)
%PARSEXML Read XML file into an XMLElement object.
    try
       tree = xmlread(filename);
    catch
       error('Failed to read XML file %s.',filename);
    end
    % Recurse over child nodes. This could run into problems
    % with very deeply nested trees.
    elt = parseChildNodes(tree);
end

% ----- Local function PARSECHILDNODES -----
function childElts = parseChildNodes(node, parentElt)
% Recurse over node children.
    if node.hasChildNodes
        childNodes = node.getChildNodes;
        numChildNodes = childNodes.getLength;
        childElts(numChildNodes) = XMLElement;
        c = 0;
        for i = 0:numChildNodes-1
            child = childNodes.item(i);
            if isa(child, 'org.apache.xerces.dom.TextImpl')
                if i == 0 && nargin == 2
                    parentElt.text = char(child.getData);
                elseif i >= 2
                    childElts(c).tail = char(child.getData);
                end
            else
                c = c + 1;
                childElts(c) = makeElementFromNode(child);
            end
        end
        childElts = childElts(1:c);
    else
        childElts = [];
    end
end

% ----- Local function MAKEELEMENTFROMNODE -----
function elt = makeElementFromNode(node)
% Create Element representation of node info.
    elt = XMLElement;
    elt.tag = char(node.getNodeName);
    elt.attributes = parseAttributes(node);
    elt.children = parseChildNodes(node, elt);
end

% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(node)
% Create attributes structure.
    attributes = [];
    if node.hasAttributes
       theAttributes = node.getAttributes;
       numAttributes = theAttributes.getLength;
       allocCell = cell(1, numAttributes);
       attributes = struct('name', allocCell, 'value', allocCell);
       for i = 1:numAttributes
          attrib = theAttributes.item(i-1);
          attributes(i).name = char(attrib.getName);
          attributes(i).value = char(attrib.getValue);
       end
    end
end

% ======================================================
