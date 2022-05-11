/**
 * @file avltree.cpp
 */
using namespace std;

template <class K, class V>
V AVLTree<K, V>::find(const K& key) const
{
    return find(root, key);
}

template <class K, class V>
V AVLTree<K, V>::find(Node* subtree, const K& key) const
{
    if (subtree == NULL)
        return V();
    else if (key == subtree->key)
        return subtree->value;
    else {
        if (key < subtree->key)
            return find(subtree->left, key);
        else
            return find(subtree->right, key);
    }
}

template <class K, class V>
void AVLTree<K, V>::rotateLeft(Node*& t)
{
    functionCalls.push_back("rotateLeft"); // Stores the rotation name (don't remove this)
    Node* a = t;
    Node* b = t->right;

    a->right = b->left;
    b->left = a;
    t= b;

    a->height = 1 + max(heightOrNeg1(a->left), heightOrNeg1(a->right));
    b->height = 1 + max(heightOrNeg1(b->left), heightOrNeg1(b->right));

}

template <class K, class V>
void AVLTree<K, V>::rotateLeftRight(Node*& t)
{
    functionCalls.push_back("rotateLeftRight"); // Stores the rotation name (don't remove this)
    rotateLeft(t->left);
    rotateRight(t);
}

template <class K, class V>
void AVLTree<K, V>::rotateRight(Node*& t)
{
    functionCalls.push_back("rotateRight"); // Stores the rotation name (don't remove this)
    Node* a = t;
    Node* b = t->left;

    a->left = b->right;
    b->right = a;
    t= b;

    a->height = 1 + max(heightOrNeg1(a->left), heightOrNeg1(a->right));
    b->height = 1 + max(heightOrNeg1(b->left), heightOrNeg1(b->right));

}

template <class K, class V>
void AVLTree<K, V>::rotateRightLeft(Node*& t)
{
    functionCalls.push_back("rotateRightLeft"); // Stores the rotation name (don't remove this)
    rotateRight(t->right);
    rotateLeft(t);
}

template <class K, class V>
void AVLTree<K, V>::rebalance(Node*& subtree)
{
    if(subtree == NULL) return;
    if (heightOrNeg1(subtree->right) - heightOrNeg1(subtree->left) < -1) {
        int leftbalance_ = heightOrNeg1(subtree->left->right) - heightOrNeg1(subtree->left->left);
        if (leftbalance_ == -1) rotateRight(subtree);
        else rotateLeftRight(subtree);
    } else if (heightOrNeg1(subtree->right) - heightOrNeg1(subtree->left) > 1) {
        int rightbalance_ = heightOrNeg1(subtree->right->right) - heightOrNeg1(subtree->right->left);
        if (rightbalance_ == 1) rotateLeft(subtree);
        else rotateRightLeft(subtree);
    }
    subtree->height =  1 + max(heightOrNeg1(subtree->right), heightOrNeg1(subtree->left));
}

template <class K, class V>
void AVLTree<K, V>::insert(const K & key, const V & value)
{
    insert(root, key, value);
}

template <class K, class V>
void AVLTree<K, V>::insert(Node*& subtree, const K& key, const V& value)
{
    if (!subtree) subtree = new Node(key, value);
    else if (key < subtree->key) insert(subtree->left, key, value);
    else if (key > subtree->key) insert(subtree->right, key, value);
    rebalance(subtree);
}

template <class K, class V>
void AVLTree<K, V>::remove(const K& key)
{
    remove(root, key);
}

template <class K, class V>
void AVLTree<K, V>::remove(Node*& subtree, const K& key)
{
    if (subtree == NULL)
        return;

    if (key < subtree->key) {
        
        remove(subtree->left, key);
    } else if (key > subtree->key) { 
        
        remove(subtree->right, key);
    } else {
        if (subtree->left == NULL && subtree->right == NULL) {
            /* no-child remove */
            
            delete subtree;
            subtree = NULL;
        } else if (subtree->left != NULL && subtree->right != NULL) {
            /* two-child remove */
            
            Node* temp = subtree->left;
            while (temp->right != NULL) {
                temp = temp->right;
            }
            swap(subtree, temp);
            remove(subtree->left, key);
            
        } else {
            /* one-child remove */
            
            Node* leftTree = subtree->left;
            Node* rightTree = subtree->right;
            if (leftTree != NULL) {
                delete subtree;
                subtree = leftTree;
            }
            if (rightTree != NULL) {
                delete subtree;
                subtree = rightTree;
            }
        }
    }
    rebalance(subtree);
}
template<class K, class V>
void AVLTree<K, V>::BFS_Traversal(AVLTree<int, string> tree, string inputFile, string outputFile)
{
    ifstream inputBFS (inputFile);
    ofstream outputBFS;
    string strBFS;
    int i = 0;
    outputBFS.open(outputFile);
    while (getline(inputBFS, strBFS)) {
        tree.insert(i, strBFS);
        i++;
    }
    vector<int> traversal = tree.getLevelorderTraversal();
    for (unsigned i = 0; i < traversal.size(); i++) {
        outputBFS << tree.find(traversal[i]) << endl;
    }
}