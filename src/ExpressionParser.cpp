#include "ExpressionParser.h"
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <iostream>

Node::Node(char v) : value(v), left(nullptr), right(nullptr) {}

Node::~Node() {
    delete left;
    delete right;
}

ExpressionParser::ExpressionParser(const std::string& input) : input(input), currentIndex(0) {
    tokenize(this->input);
}

std::vector<Token> ExpressionParser::tokenize(const std::string& input) {
    tokens.clear();
    int open_pars = 0;
    for (char c : input) {
        switch (c) {
        case 'U':
            tokens.push_back(Token(UNION, c));
            break;
        case 'I':
            tokens.push_back(Token(INTERSECT, c));
            break;
        case '\\':
            tokens.push_back(Token(DIFFERENCE, c));
            break;
        case '(':
            tokens.push_back(Token(OPEN_PAR, c));
            open_pars++;
            break;
        case ')':
            tokens.push_back(Token(CLOSE_PAR, c));
            open_pars--;
            break;
        case 'N':
            tokens.push_back(Token(NAND, c));
            break;
        case 'X':
            tokens.push_back(Token(XOR, c));
            break;
        default:
            tokens.push_back(Token(VALUE, c));
            break;
        }
    }
    if (open_pars != 0) {
        throw std::runtime_error("[WARNING]: Unmatched parentheses in input.");
    }
    return tokens;
}

Node* ExpressionParser::parse() {
    currentIndex = 0;
    Node* result = parseExpression();
    if (currentIndex != tokens.size()) {
        throw std::runtime_error("[WARNING]: Unexpected tokens after parsing.");
    }
    return result;
}

Node* ExpressionParser::parseExpression() {
    return parseUnionIntersection();
}

Node* ExpressionParser::parseUnionIntersection() {
    Node* node = parseDifferenceNandXor();
    while (currentIndex < tokens.size() &&
        (tokens[currentIndex].type == UNION || tokens[currentIndex].type == INTERSECT)) {
        Token op = tokens[currentIndex++];
        Node* right = parseDifferenceNandXor();
        node = createNodeForToken(op, node, right);
    }
    return node;
}

Node* ExpressionParser::parseDifferenceNandXor() {
    Node* node = parsePrimary();
    while (currentIndex < tokens.size() &&
        (tokens[currentIndex].type == DIFFERENCE || tokens[currentIndex].type == NAND || tokens[currentIndex].type == XOR)) {
        Token op = tokens[currentIndex++];
        Node* right = parsePrimary();
        node = createNodeForToken(op, node, right);
    }
    return node;
}

Node* ExpressionParser::parsePrimary() {
    if (currentIndex >= tokens.size()) throw std::runtime_error("[WARNING]: Incomplete expression.");
    Token token = tokens[currentIndex];

    if (token.type == OPEN_PAR) {
        currentIndex++; // Skip '('
        Node* node = parseExpression();
        if (currentIndex >= tokens.size() || tokens[currentIndex].type != CLOSE_PAR) {
            throw std::runtime_error("[WARNING]: Missing closing parenthesis.");
        }
        currentIndex++; // Skip ')'
        return node;
    }
    else if (token.type == VALUE) {
        currentIndex++;
        return new Node(token.value);
    }
    else {
        throw std::runtime_error("[WARNING]: Unexpected token.");
    }
}

Node* ExpressionParser::createNodeForToken(const Token& token, Node* left, Node* right) {
    Node* node = new Node(token.value);
    node->left = left;
    node->right = right;
    return node;
}