#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

#include <vector>
#include <string>
#include <map>
#include <set>
#include <sstream>

enum TokenType {
    VALUE,
    UNION,
    INTERSECT,
    DIFFERENCE,
    NAND,
    XOR,
    OPEN_PAR,
    CLOSE_PAR
};

struct Token {
    TokenType type;
    char value;
    Token(TokenType t, char v) : type(t), value(v) {}
};

struct Node {
    char value;
    Node* left;
    Node* right;

    Node(char v);
    ~Node();
};

class ExpressionParser {
public:
    ExpressionParser(const std::string& input);
    Node* parse();
    std::vector<Token> tokenize(const std::string& input);

private:
    std::string input;
    std::vector<Token> tokens;
    size_t currentIndex;

    Node* parseExpression();
    Node* parseUnionIntersection();
    Node* parseDifferenceNandXor();
    Node* parsePrimary();
    Node* createNodeForToken(const Token& token, Node* left, Node* right);
};

#endif // EXPRESSION_PARSER_H
