import java.util.*;
// import java.util.Stack;

class HashMap {

    // 895. Maximum Frequency Stack
    /*
     * Time Complexity : O(n) Space Complexity: O(n)
     * 
     * Use two hashmaps -> one for stack and one for freq
     * In stack Map we maintain
     * {freq : stack} So, whenever an element is pushed in stack, we check its
     * frequency from freq map, then in stack map, we push it in the stack of that
     * freq
     */
    class FreqStack {
        HashMap<Integer, Stack<Integer>> st;
        HashMap<Integer, Integer> freqMap;
        int maxFreq;

        public FreqStack() {
            st = new HashMap<>();
            freqMap = new HashMap<>();
            maxFreq = 0;
        }

        public void push(int val) {
            // get the freq of val and increase by 1
            int freq = freqMap.getOrDefault(val, 0) + 1;

            // update freq in freqMap
            freqMap.put(val, freq);

            // push in stack in stack map
            if (!st.containsKey(freq))
                st.put(freq, new Stack<Integer>());
            st.get(freq).push(val);

            // update maxFreq
            maxFreq = Math.max(maxFreq, freq);
        }

        public int pop() {
            // pop element from stack map
            int res = st.get(maxFreq).pop();

            if (st.get(maxFreq).size() == 0) {
                st.remove(maxFreq);
                maxFreq--;
            }

            // update freq map
            int freq = freqMap.get(res) - 1;
            freqMap.put(res, freq);

            if (freq == 0)
                freqMap.remove(res);

            return res;
        }
    }

    /**
     * Your FreqStack object will be instantiated and called as such: FreqStack obj
     * = new FreqStack(); obj.push(val); int param_2 = obj.pop();
     */
    public static void main(String[] args) {

    }
}