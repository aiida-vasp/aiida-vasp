"""
Wrapper class to simply updating an Node when configuring inputs of a calculation.
"""
from collections import UserDict
from typing import Any

from aiida.orm.nodes.data.dict import Dict


class DictWrapper(UserDict):
    """
    Wrapper class for a dict

    Wraps around an Dict node, saved or not, to allow python dict-like in place modifications.
    """

    def __init__(self, node, namespace=None, port=None):
        """
        Wrapper around a ``Dict`` node. Optionally also include a namespace and a port that the
        underlying Node should be assigned to.
        The wrapped ``Dict`` node is never update in practice and only used as an initial reference data.
        If there is no change, the same Dict node will be passed as it is under the ``node`` property.

        :param node: The node to be wrapped.
        :param namespace: The namespace that the node should be assigned to.
        :param port: The port under the namespace.

        A new node is created whenever necessary.
        """
        super().__init__()

        self._stored_node = node
        self._unstored_node = None  # A new node to be used
        self.namespace = namespace
        self.port = port
        self.data = node.get_dict()
        if namespace is not None:
            namespace[port] = node

    def __getattr__(self, name: str) -> Any:
        return getattr(self.node, name)

    @property
    def is_updated(self):
        return self.data != self._stored_node.get_dict()

    @property
    def node(self):
        """
        Returns a node object that represents the stored data.

        This can be a stored node (if no update) or an unstore node
        """
        if self.is_updated:
            if self._unstored_node is None:
                self._unstored_node = Dict(dict=self.data)
            return self._unstored_node
        return self._stored_node

    def __setitem__(self, key, value):
        """
        Set the value of a key
        """
        update = False
        if key in self.data:
            if self.data[key] != value:
                update = True
        else:
            update = True
        self._ensure_unstored()
        super().__setitem__(key, value)
        # If we have updated the value of a key
        if update:
            if self._unstored_node is None:
                self._unstored_node = Dict(dict=self.data)
            else:
                self._unstored_node.base.attributes.set(key, value)

        if self.namespace is not None:
            self.namespace[self.port] = self.node

    def _ensure_unstored(self):
        """
        Ensure that self._unstored_node is indeed unstored
        If it is not the case - we move it to self._stored_node an create a new
        unstored node.
        """
        if self._unstored_node and self._unstored_node.is_stored:
            self._stored_node = self._unstored_node
            self._unstored_node = Dict(dict=self.data)

    def __delitem__(self, key) -> None:
        """
        Delete an item

        Delete from both the python dictionary and the underlying node
        """

        self._ensure_unstored()
        super().__delitem__(key)
        if self._unstored_node is None:
            self._unstored_node = Dict(dict=self.data)
        else:
            self._unstored_node.delete_attribute(key)

        if self.namespace is not None:
            self.namespace[self.port] = self.node

    def validate(self):
        """Validate the consistency between the node and the python dictionary"""
        assert self.data == self.node.get_dict()

    @classmethod
    def serializer(cls, data, port=None):
        """
        The `serializer` function can be used for `ProcessSpec.input`.

        Returns a Dict node for the given data.
        """
        _ = port
        if isinstance(data, DictWrapper):
            return data.node
        return Dict(dict=data)
