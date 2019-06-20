#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <QVector>
#include <QString>
#include <QObject>

using namespace std;

class codepart{ //main class
public:
    QString organization;
    QString project;
    QString component;
    QString cllass;
    int method;
    int creator;

    codepart(){
    }

    codepart(QString cnt, QString reg, QString cit, QString str, int hos, int rom){
        organization = cnt;
        project = reg;
        component = cit;
        cllass = str;
        method = hos;
        creator = rom;
    }

    friend bool operator> (const codepart &d1, const codepart &d2);
    friend bool operator<= (const codepart &d1, const codepart &d2);
    friend bool operator< (const codepart &d1, const codepart &d2);
    friend bool operator>= (const codepart &d1, const codepart &d2);
    QString getRegion() const;
    void setRegion(const QString &value);
};

bool operator> (const codepart &d1, const codepart &d2){
    return d1.creator > d2.creator;
}

bool operator>= (const codepart &d1, const codepart &d2){
    return d1.creator >= d2.creator;
}

bool operator< (const codepart &d1, const codepart &d2){
    return d1.creator < d2.creator;
}

bool operator<= (const codepart &d1, const codepart &d2){
    return d1.creator <= d2.creator;
}

class database{
public:
    QVector<codepart> data;

    database(){
    }

    void add_rand(){
        QString cnt = QString::fromStdString("Organization" + to_string(data.size()));
        QString reg = QString::fromStdString("Project" + to_string(data.size()));
        QString cit = QString::fromStdString("Component" + to_string(data.size()));
        QString str = QString::fromStdString("Class" + to_string(data.size()));
        int hos = rand()%100 + 1;
        int rom = rand()%200 + 1;;
        data.push_back(codepart(cnt, reg, cit, str, hos, rom));
    }

    void add(QString cnt, QString reg, QString cit, QString str, int hos, int rom){
        data.push_back(codepart(cnt, reg, cit, str, hos, rom));
    }
};

/* list node */
struct Node{
    codepart * a;
    Node * next;
};

class MyList{ //SL circular
public:
    Node * head;
    MyList(){
        head = nullptr;
    }
    ~MyList(){
        while (head != nullptr){
            Node *tmp = head->next;
            delete head;
            head = tmp;
        }
    }

    void add(codepart * x){ //adding
        Node * tmp = new Node;
        tmp->a = x;
        tmp->next = head;
        head = tmp;
    }

    void del(codepart * x){ //deleting
        if (head->a == x){
            Node *tmp = head->next;
            delete head;
            head = tmp;
            return;
        } else {
            for (auto i = head; i != nullptr; i = i->next){
                Node * k = i->next;
                if (k->a == x){
                    i->next = k->next;
                    delete k;
                    return;
                }
            }
        }
    }

    bool implement(codepart * x){ //implemetion
        for (auto i = head; i != nullptr; i = i->next){
            if (i->a == x){
                return true;
            }
        }
        return false;
    }

    int size(){
        int k = 0;
        for (auto i = head; i != nullptr; i = i->next){
            k++;
        }
        return k;
    }
};

// An AVL tree node
class node{
    public:
    codepart * key;
    node *left, *right;
};

/* Helper function that allocates a new node with the given key and NULL left and right pointers.*/
node* newNode(codepart * key){
    node* Node = new node();
    Node->key = key;
    Node->left = Node->right = nullptr;
    return (Node);
}

/* A utility function to right rotate subtree rooted with y.*/
node *rightRotate(node *x){
    node *y = x->left;
    x->left = y->right;
    y->right = x;
    return y;
}

/* A utility function to left rotate subtree rooted with x.*/
node *leftRotate(node *x){
    node *y = x->right;
    x->right = y->left;
    y->left = x;
    return y;
}

/* This function brings the key at root if key is present in tree. If key is not present, then it brings the last accessed item at
root. This function modifies the tree and returns the new root.*/
node *splay(node *root, codepart * key)
{
    // Base cases: root is NULL or
    // key is present at root
    if (root == nullptr || root->key == key)
        return root;

    // Key lies in left subtree
    if (*root->key > *key){
        // Key is not in tree, we are done
        if (root->left == nullptr) return root;

        // Zig-Zig (Left Left)
        if (*root->left->key > *key){
            // First recursively bring the
            // key as root of left-left
            root->left->left = splay(root->left->left, key);

            // Do first rotation for root,
            // second rotation is done after else
            root = rightRotate(root);
        }
        else if (*root->left->key < *key){ // Zig-Zag (Left Right)
            // First recursively bring
            // the key as root of left-right
            root->left->right = splay(root->left->right, key);

            // Do first rotation for root->left
            if (root->left->right != nullptr)
                root->left = leftRotate(root->left);
        }

        // Do second rotation for root
        return (root->left == nullptr)? root: rightRotate(root);
    }
    else // Key lies in right subtree
    {
        // Key is not in tree, we are done
        if (root->right == nullptr) return root;

        // Zag-Zig (Right Left)
        if (*root->right->key > *key)
        {
            // Bring the key as root of right-left
            root->right->left = splay(root->right->left, key);

            // Do first rotation for root->right
            if (root->right->left != nullptr)
                root->right = rightRotate(root->right);
        }
        else if (*root->right->key < *key)// Zag-Zag (Right Right)
        {
            // Bring the key as root of
            // right-right and do first rotation
            root->right->right = splay(root->right->right, key);
            root = leftRotate(root);
        }

        // Do second rotation for root
        return (root->right == nullptr)? root: leftRotate(root);
    }
}

// Function to insert a new key k in splay tree with given root
node *insert(node *root, codepart * k){
    // Simple Case: If tree is empty
    if (root ==nullptr) return newNode(k);

    // Bring the closest leaf node to root
    root = splay(root, k);

    // If key is already present, then return
    if (root->key == k) return root;

    // Otherwise allocate memory for new node
    node *newnode = newNode(k);

    // If root's key is greater, make root as right child of newnode and copy the left child of root to newnode
    if (*root->key > *k){
        newnode->right = root;
        newnode->left = root->left;
        root->left = nullptr;
    }
    // If root's key is smaller, make root as left child of newnode and copy the right child of root to newnode
    else {
        newnode->left = root;
        newnode->right = root->right;
        root->right = nullptr;
    }

    return newnode; // newnode becomes new root
}

// The delete function for Splay tree. Note that this function returns the new root of Splay Tree after removing the key
node* delete_key(node * root, codepart * key){
    node *temp;
    if (!root)
        return nullptr;

    // Splay the given key
    root = splay(root, key);

    // If key is not present, then
    // return root
    if (key != root->key)
        return root;

    // If key is present
    // If left child of root does not exist
    // make root->right as root
    if (!root->left){
        temp = root;
        root = root->right;
    }
    // Else if left child exits
    else {
        temp = root;
        /*Note: Since key == root->key,
        so after Splay(key, root->lchild),
        the tree we get will have no right child tree
        and maximum node in left subtree will get splayed*/
        // New root
        root = splay(root->left, key);

        // Make right child of previous root  as
        // new root's right child
        root->right = temp->right;
    }

    // free the previous root node, that is,
    // the node containing the key
    free(temp);

    // return root of the new Splay Tree
    return root;
}

/* The search function for Splay tree. Note that this function returns the new root of Splay Tree. If key is present in tree then, it is moved to root.*/
node *search(node *root, codepart * key)
{
    return splay(root, key);
}

//BTree
class BTreeNode
{
    int *keys;  // An array of keys
    int t;      // Minimum degree (defines the range for number of keys)
    BTreeNode **C; // An array of child pointers
    int n;     // Current number of keys
    bool leaf; // Is true when node is leaf. Otherwise false
public:
    BTreeNode(int _t, bool _leaf);   // Constructor

    // A utility function to insert a new key in the subtree rooted with
    // this node. The assumption is, the node must be non-full when this
    // function is called
    void insertNonFull(int k);

    // A utility function to split the child y of this node. i is index of y in
    // child array C[].  The Child y must be full when this function is called
    void splitChild(int i, BTreeNode *y);

    // A function to traverse all nodes in a subtree rooted with this node
    void traverse();

    // A function to search a key in subtree rooted with this node.
    BTreeNode *search(int k);   // returns NULL if k is not present.

// Make BTree friend of this so that we can access private members of this
// class in BTree functions
friend class BTree;
};

// A BTree
class BTree
{
    BTreeNode *root; // Pointer to root node
    int t;  // Minimum degree
public:
    // Constructor (Initializes tree as empty)
    BTree(int _t)
    {  root = NULL;  t = _t; }

    // function to traverse the tree
    void traverse()
    {  if (root != NULL) root->traverse(); }

    // function to search a key in this tree
    BTreeNode* search(int k)
    {  return (root == NULL)? NULL : root->search(k); }

    // The main function that inserts a new key in this B-Tree
    void insert(int k);
};

// Constructor for BTreeNode class
BTreeNode::BTreeNode(int t1, bool leaf1)
{
    // Copy the given minimum degree and leaf property
    t = t1;
    leaf = leaf1;

    // Allocate memory for maximum number of possible keys
    // and child pointers
    keys = new int[2*t-1];
    C = new BTreeNode *[2*t];

    // Initialize the number of keys as 0
    n = 0;
}

// Function to search key k in subtree rooted with this node
BTreeNode *BTreeNode::search(int k)
{
    // Find the first key greater than or equal to k
    int i = 0;
    while (i < n && k > keys[i])
        i++;

    // If the found key is equal to k, return this node
    if (keys[i] == k)
        return this;

    // If key is not found here and this is a leaf node
    if (leaf == true)
        return NULL;

    // Go to the appropriate child
    return C[i]->search(k);
}

// The main function that inserts a new key in this B-Tree
void BTree::insert(int k)
{
    // If tree is empty
    if (root == NULL)
    {
        // Allocate memory for root
        root = new BTreeNode(t, true);
        root->keys[0] = k;  // Insert key
        root->n = 1;  // Update number of keys in root
    }
    else // If tree is not empty
    {
        // If root is full, then tree grows in height
        if (root->n == 2*t-1)
        {
            // Allocate memory for new root
            BTreeNode *s = new BTreeNode(t, false);

            // Make old root as child of new root
            s->C[0] = root;

            // Split the old root and move 1 key to the new root
            s->splitChild(0, root);

            // New root has two children now.  Decide which of the
            // two children is going to have new key
            int i = 0;
            if (s->keys[0] < k)
                i++;
            s->C[i]->insertNonFull(k);

            // Change root
            root = s;
        }
        else  // If root is not full, call insertNonFull for root
            root->insertNonFull(k);
    }
}

// A utility function to insert a new key in this node
// The assumption is, the node must be non-full when this
// function is called
void BTreeNode::insertNonFull(int k)
{
    // Initialize index as index of rightmost element
    int i = n-1;

    // If this is a leaf node
    if (leaf == true)
    {
        // The following loop does two things
        // a) Finds the location of new key to be inserted
        // b) Moves all greater keys to one place ahead
        while (i >= 0 && keys[i] > k)
        {
            keys[i+1] = keys[i];
            i--;
        }

        // Insert the new key at found location
        keys[i+1] = k;
        n = n+1;
    }
    else // If this node is not leaf
    {
        // Find the child which is going to have the new key
        while (i >= 0 && keys[i] > k)
            i--;

        // See if the found child is full
        if (C[i+1]->n == 2*t-1)
        {
            // If the child is full, then split it
            splitChild(i+1, C[i+1]);

            // After split, the middle key of C[i] goes up and
            // C[i] is splitted into two.  See which of the two
            // is going to have the new key
            if (keys[i+1] < k)
                i++;
        }
        C[i+1]->insertNonFull(k);
    }
}

// A utility function to split the child y of this node
// Note that y must be full when this function is called
void BTreeNode::splitChild(int i, BTreeNode *y)
{
    // Create a new node which is going to store (t-1) keys
    // of y
    BTreeNode *z = new BTreeNode(y->t, y->leaf);
    z->n = t - 1;

    // Copy the last (t-1) keys of y to z
    for (int j = 0; j < t-1; j++)
        z->keys[j] = y->keys[j+t];

    // Copy the last t children of y to z
    if (y->leaf == false)
    {
        for (int j = 0; j < t; j++)
            z->C[j] = y->C[j+t];
    }

    // Reduce the number of keys in y
    y->n = t - 1;

    // Since this node is going to have a new child,
    // create space of new child
    for (int j = n; j >= i+1; j--)
        C[j+1] = C[j];

    // Link the new child to this node
    C[i+1] = z;

    // A key of y will move to this node. Find location of
    // new key and move all greater keys one space ahead
    for (int j = n-1; j >= i; j--)
        keys[j+1] = keys[j];

    // Copy the middle key of y to this node
    keys[i] = y->keys[t-1];

    // Increment count of keys in this node
    n = n + 1;
}

//Hash Table with Linear Probing
//template for generic type
template<typename K, typename V>

//Hashnode class
class hashnode{
    public:
    V value;
    K key;

    //Constructor of hashnode
    hashnode(K key, V value)
    {
        this->value = value;
        this->key = key;
    }
};

//template for generic type
template<typename K, typename V>

//Our own Hashmap class
class hashmap
{
public:
    //hash element array
    hashnode<K,V> **arr;
    int capacity;
    //current size
    int size;
    //dummy node
    hashnode<K,V> *dummy;
    hashmap()
    {
        //Initial capacity of hash array
        capacity = 128;
        size=0;
        arr = new hashnode<K,V>*[capacity];

        //Initialise all elements of array as NULL
        for(int i=0 ; i < capacity ; i++)
            arr[i] = NULL;

        //dummy node with value and key -1
        dummy = new hashnode<K,V>(-1, nullptr);
    }
    // This implements hash function to find index
    // for a key
    int hashCode(K key)
    {
        return key % capacity;
    }

    //Function to add key value pair
    void insertNode(K key, V value)
    {
        hashnode<K,V> *temp = new hashnode<K,V>(key, value);

        // Apply hash function to find index for given key
        int hashIndex = hashCode(key);

        //find next free space
        while(arr[hashIndex] != NULL && arr[hashIndex]->key != key
               && arr[hashIndex]->key != -1)
        {
            hashIndex++;
            hashIndex %= capacity;
        }

        //if new node to be inserted increase the current size
        if(arr[hashIndex] == NULL || arr[hashIndex]->key == -1)
            size++;
        arr[hashIndex] = temp;
    }

    //Function to delete a key value pair
    V deleteNode(int key)
    {
        // Apply hash function to find index for given key
        int hashIndex = hashCode(key);

        //finding the node with given key
        while(arr[hashIndex] != NULL)
        {
            //if node found
            if(arr[hashIndex]->key == key)
            {
                hashnode<K,V> *temp = arr[hashIndex];

                //Insert dummy node here for further use
                arr[hashIndex] = dummy;

                // Reduce size
                size--;
                return temp->value;
            }
            hashIndex++;
            hashIndex %= capacity;

        }

        //If not found return null
        return NULL;
    }

    //Function to search the value for a given key
    V get(int key)
    {
        // Apply hash function to find index for given key
        int hashIndex = hashCode(key);
        //finding the node with given key
        while(arr[hashIndex] != NULL)
        {    int counter =0;
             if(counter++>capacity)  //to avoid infinite loop
                return NULL;
            //if node found return its value
            if(arr[hashIndex]->key == key)
                return arr[hashIndex]->value;
            hashIndex++;
            hashIndex %= capacity;
        }

        //If not found return null
        return NULL;
    }

    //Return current size
    int sizeofMap()
    {
        return size;
    }

    //Return true if size is 0
    bool isEmpty()
    {
        return size == 0;
    }
};

//Hash Table chaining
// HashNode Class Declaration
const int TABLE_SIZE = 128;

class HashNode
{
    public:
    int key;
    codepart * value;
    HashNode* next;
    HashNode(int key, codepart * value){
        this->key = key;
        this->value = value;
        this->next = nullptr;
    }
};

//HashMap Class Declaration
class HashMap
{
    private:
        HashNode** htable;
    public:
        HashMap()
        {
            htable = new HashNode*[TABLE_SIZE];
            for (int i = 0; i < TABLE_SIZE; i++)
                htable[i] = nullptr;
        }
        ~HashMap()
        {
            delete[] htable;
        }
        //Hash Function
        int HashFunc(int key)
        {
            return key % TABLE_SIZE;
        }

        //Insert Element at a key
        void Insert(int key, codepart * value)
        {
            int hash_val = HashFunc(key);
            HashNode* prev = nullptr;
            HashNode* entry = htable[hash_val];
            while (entry != nullptr)
            {
                prev = entry;
                entry = entry->next;
            }
            if (entry == nullptr)
            {
                entry = new HashNode(key, value);
                if (prev == nullptr)
            {
                    htable[hash_val] = entry;
                }
            else
            {
                    prev->next = entry;
                }
            }
            else
            {
                entry->value = value;
            }
        }

        //Remove Element at a key
        void Remove(int key)
        {
            int hash_val = HashFunc(key);
            HashNode* entry = htable[hash_val];
            HashNode* prev = nullptr;
            if (entry == nullptr || entry->key != key)
            {
                return;
            }
            while (entry->next != nullptr)
        {
                prev = entry;
                entry = entry->next;
            }
            if (prev != nullptr)
            {
                prev->next = entry->next;
            }
            delete entry;
        }
        //Search Element at a key
        codepart * Search(int key)
        {
            int hash_val = HashFunc(key);
            HashNode* entry = htable[hash_val];
            while (entry != nullptr){
                if (entry->key == key){
                    return entry->value;
                }
                entry = entry->next;
            }
                return nullptr;
        }
};

//insertion sort
void insertionSort(QVector<int> &vec){
    for (auto it = vec.begin(); it != vec.end(); it++)
    {
        // Searching the upper bound, i.e., first
        // element greater than *it from beginning
        auto const insertion_point = upper_bound(vec.begin(), it, *it);

        // Shifting the unsorted part
        rotate(insertion_point, it, it+1);
    }
}

//heapsort
void Swap(QVector<int>& vHeap, QVector<int>::size_type i, QVector<int>::size_type j)
{
    if(i == j)
        return;

    int temp;
    temp = vHeap[i];
    vHeap[i] = vHeap[j];
    vHeap[j] = temp;
}

void Sift(QVector<int>& vHeap, const QVector<int>::size_type heapSize, const QVector<int>::size_type siftNode)
{
    QVector<int>::size_type i, j;

    j = siftNode;
    do
    {
        i = j;
        if(((2*i + 1) < heapSize) && vHeap[j] < vHeap[2*i + 1])
            j = 2*i + 1;
        if(((2*i + 2) < heapSize) && vHeap[j] < vHeap[2*i + 2])
            j = 2*i + 2;

        Swap(vHeap, i, j);
    }
    while(i != j);
}

void MakeInitialHeap(QVector<int>& vHeap)
{
    for(int i = vHeap.size() - 1; i >= 0; --i)
    {
        Sift(vHeap, vHeap.size(), i);
    }
}

void HeapSort(QVector<int>& vHeap)
{
    MakeInitialHeap(vHeap);
    for(std::vector<int>::size_type i = vHeap.size()-1; i > 0; --i)
    {
        Swap(vHeap, i, 0);
        Sift(vHeap, i, 0);
    }
}

//merge sort
QVector<int> merge(QVector<int> left, QVector<int> right)
{
    QVector<int> result;
    while (left.size() > 0 || right.size() > 0) {
        if (left.size() > 0 && right.size() > 0) {
            if (left.front() <= right.front()) {
                result.push_back(left.front());
                left.erase(left.begin());
            }
            else {
                result.push_back(right.front());
                right.erase(right.begin());
            }
        }  else if (left.size() > 0) {
            for (int i = 0; i < left.size(); i++)
                result.push_back(left[i]);
            break;
        }  else if (right.size() > 0) {
            for (int i = 0; i < right.size(); i++)
                result.push_back(right[i]);
            break;
        }
    }
    return result;
}

QVector<int> mSort(QVector<int> m){
    if (m.size() <= 1)
        return m;

    QVector<int> left, right, result;
    int middle = (m.size()+ 1) / 2;

    for (int i = 0; i < middle; i++) {
        left.push_back(m[i]);
    }

    for (int i = middle; i < m.size(); i++) {
        right.push_back(m[i]);
    }

    left = mSort(left);
    right = mSort(right);
    result = merge(left, right);

    return result;
}

void MergeSort(QVector<int> &m){
    m = mSort(m);
}

//count sort
void countSort(QVector<int> &a) {
    int d = 8, w = 32;
    for (int p = 0; p < w/d; p++) {
        QVector<int> c(1<<d, 0);
        // the next three for loops implement counting-sort
        QVector<int> b(a.size());
        for (int i = 0; i < a.size(); i++)
            c[(a[i] >> d*p)&((1<<d)-1)]++;
        for (int i = 1; i < 1<<d; i++)
            c[i] += c[i-1];
        for (int i = a.size()-1; i >= 0; i--)
            b[--c[(a[i] >> d*p)&((1<<d)-1)]] = a[i];
        a = b;
    }
}

#endif // ALGORITHMS_H

QString codepart::getRegion() const
{
return project;
}

void codepart::setRegion(const QString &value)
{
project = value;
}
